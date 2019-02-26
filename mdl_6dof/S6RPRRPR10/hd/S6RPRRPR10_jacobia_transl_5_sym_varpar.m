% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR10_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR10_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:03
% EndTime: 2019-02-26 21:06:03
% DurationCPUTime: 0.10s
% Computational Cost: add. (64->26), mult. (144->33), div. (0->0), fcn. (167->6), ass. (0->23)
t10 = cos(qJ(3));
t17 = pkin(8) + r_i_i_C(2);
t14 = r_i_i_C(3) + qJ(5);
t18 = pkin(4) + r_i_i_C(1);
t6 = sin(qJ(4));
t9 = cos(qJ(4));
t22 = t14 * t6 + t18 * t9 + pkin(3);
t7 = sin(qJ(3));
t25 = t22 * t10 + t17 * t7;
t23 = t17 * t10;
t8 = sin(qJ(1));
t21 = t8 * t6;
t20 = t8 * t9;
t19 = pkin(1) + pkin(7);
t11 = cos(qJ(1));
t16 = t11 * t6;
t15 = t11 * t9;
t13 = t7 * pkin(3) + qJ(2) - t23;
t4 = t7 * t15 - t21;
t3 = t7 * t16 + t20;
t2 = t7 * t20 + t16;
t1 = t7 * t21 - t15;
t5 = [t13 * t11 + t14 * t3 + t18 * t4 - t19 * t8, t8, t25 * t8, -t18 * t1 + t14 * t2, t1, 0; t14 * t1 + t19 * t11 + t13 * t8 + t18 * t2, -t11, -t25 * t11, -t14 * t4 + t18 * t3, -t3, 0; 0, 0, -t22 * t7 + t23 (t14 * t9 - t18 * t6) * t10, t10 * t6, 0;];
Ja_transl  = t5;
