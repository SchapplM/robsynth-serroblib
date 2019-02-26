% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP8
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
% Datum: 2019-02-26 21:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPP8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_jacobia_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:00:20
% EndTime: 2019-02-26 21:00:20
% DurationCPUTime: 0.11s
% Computational Cost: add. (83->27), mult. (186->34), div. (0->0), fcn. (218->6), ass. (0->23)
t11 = cos(qJ(3));
t16 = pkin(5) + pkin(8) + r_i_i_C(1);
t10 = cos(qJ(4));
t15 = pkin(4) + r_i_i_C(3) + qJ(6);
t17 = r_i_i_C(2) + qJ(5);
t7 = sin(qJ(4));
t23 = t15 * t10 + t17 * t7 + pkin(3);
t8 = sin(qJ(3));
t26 = t23 * t11 + t16 * t8;
t24 = t16 * t11;
t9 = sin(qJ(1));
t22 = t9 * t7;
t21 = pkin(1) + pkin(7);
t12 = cos(qJ(1));
t20 = t12 * t7;
t19 = t9 * t10;
t18 = t12 * t10;
t14 = pkin(3) * t8 + qJ(2) - t24;
t4 = t8 * t18 - t22;
t3 = t8 * t20 + t19;
t2 = t8 * t19 + t20;
t1 = t8 * t22 - t18;
t5 = [t14 * t12 + t15 * t4 + t17 * t3 - t21 * t9, t9, t26 * t9, -t15 * t1 + t17 * t2, t1, t2; t17 * t1 + t21 * t12 + t14 * t9 + t15 * t2, -t12, -t26 * t12, t15 * t3 - t17 * t4, -t3, -t4; 0, 0, -t23 * t8 + t24 (t17 * t10 - t15 * t7) * t11, t11 * t7, t11 * t10;];
Ja_transl  = t5;
