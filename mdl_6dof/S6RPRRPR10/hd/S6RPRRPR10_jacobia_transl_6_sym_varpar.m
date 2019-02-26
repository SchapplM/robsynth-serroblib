% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_transl = S6RPRRPR10_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR10_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:03
% EndTime: 2019-02-26 21:06:03
% DurationCPUTime: 0.18s
% Computational Cost: add. (122->43), mult. (293->60), div. (0->0), fcn. (358->8), ass. (0->32)
t14 = sin(qJ(3));
t18 = cos(qJ(3));
t27 = r_i_i_C(3) + pkin(9) - pkin(8);
t13 = sin(qJ(4));
t17 = cos(qJ(4));
t16 = cos(qJ(6));
t12 = sin(qJ(6));
t33 = pkin(4) + pkin(5);
t25 = -t12 * r_i_i_C(2) + t33;
t22 = t16 * r_i_i_C(1) + t25;
t26 = t12 * r_i_i_C(1) + qJ(5);
t23 = t16 * r_i_i_C(2) + t26;
t34 = t23 * t13 + t22 * t17 + pkin(3);
t37 = -t27 * t14 + t34 * t18;
t35 = t27 * t18;
t32 = -pkin(7) - pkin(1);
t19 = cos(qJ(1));
t29 = t19 * t13;
t15 = sin(qJ(1));
t30 = t15 * t17;
t8 = t14 * t29 + t30;
t3 = t8 * t16;
t31 = t15 * t13;
t28 = t19 * t17;
t24 = (-(t12 * t17 - t13 * t16) * r_i_i_C(1) - (t12 * t13 + t16 * t17) * r_i_i_C(2)) * t18;
t21 = t14 * pkin(3) + qJ(2) + t35;
t9 = t14 * t28 - t31;
t7 = t14 * t30 + t29;
t6 = t14 * t31 - t28;
t2 = t6 * t12 + t7 * t16;
t1 = -t7 * t12 + t6 * t16;
t4 = [t3 * r_i_i_C(2) + t32 * t15 + t21 * t19 + t22 * t9 + t26 * t8, t15, t37 * t15, -t22 * t6 + t23 * t7, t6, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t6 * qJ(5) + t21 * t15 - t32 * t19 + t33 * t7, -t19, -t37 * t19, t3 * r_i_i_C(1) - t23 * t9 + t25 * t8, -t8 (t9 * t12 - t3) * r_i_i_C(1) + (t8 * t12 + t9 * t16) * r_i_i_C(2); 0, 0, -t14 * t34 - t35 (qJ(5) * t17 - t33 * t13) * t18 - t24, t18 * t13, t24;];
Ja_transl  = t4;
