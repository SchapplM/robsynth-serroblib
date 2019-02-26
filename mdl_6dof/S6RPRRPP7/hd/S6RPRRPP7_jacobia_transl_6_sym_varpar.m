% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPP7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_jacobia_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:59:51
% EndTime: 2019-02-26 20:59:51
% DurationCPUTime: 0.11s
% Computational Cost: add. (83->28), mult. (181->35), div. (0->0), fcn. (211->6), ass. (0->23)
t10 = cos(qJ(3));
t14 = pkin(8) - r_i_i_C(3) - qJ(6);
t15 = pkin(4) + pkin(5) + r_i_i_C(1);
t16 = r_i_i_C(2) + qJ(5);
t6 = sin(qJ(4));
t9 = cos(qJ(4));
t22 = t15 * t9 + t16 * t6 + pkin(3);
t7 = sin(qJ(3));
t25 = t22 * t10 + t14 * t7;
t23 = t14 * t10;
t8 = sin(qJ(1));
t21 = t8 * t6;
t20 = t8 * t9;
t19 = pkin(1) + pkin(7);
t11 = cos(qJ(1));
t18 = t11 * t6;
t17 = t11 * t9;
t13 = pkin(3) * t7 + qJ(2) - t23;
t4 = t7 * t17 - t21;
t3 = t7 * t18 + t20;
t2 = t7 * t20 + t18;
t1 = t7 * t21 - t17;
t5 = [t13 * t11 + t15 * t4 + t16 * t3 - t19 * t8, t8, t25 * t8, -t15 * t1 + t16 * t2, t1, t8 * t10; t16 * t1 + t19 * t11 + t13 * t8 + t15 * t2, -t11, -t25 * t11, t15 * t3 - t16 * t4, -t3, -t11 * t10; 0, 0, -t22 * t7 + t23 (-t15 * t6 + t16 * t9) * t10, t10 * t6, -t7;];
Ja_transl  = t5;
