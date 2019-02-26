% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:01
% EndTime: 2019-02-26 22:07:01
% DurationCPUTime: 0.18s
% Computational Cost: add. (190->44), mult. (335->65), div. (0->0), fcn. (405->10), ass. (0->32)
t22 = cos(qJ(2));
t19 = sin(qJ(2));
t33 = -r_i_i_C(3) - pkin(9) - qJ(5) + pkin(8);
t30 = t33 * t19;
t39 = t22 * pkin(2) + pkin(1) + t30;
t18 = sin(qJ(3));
t21 = cos(qJ(3));
t15 = pkin(10) + qJ(6);
t14 = cos(t15);
t13 = sin(t15);
t31 = sin(pkin(10)) * pkin(5) + qJ(4);
t27 = -t13 * r_i_i_C(1) - t31;
t25 = t14 * r_i_i_C(2) - t27;
t37 = pkin(3) + cos(pkin(10)) * pkin(5) + pkin(4);
t29 = t13 * r_i_i_C(2) - t37;
t26 = t14 * r_i_i_C(1) - t29;
t38 = t18 * t25 + t21 * t26 + pkin(2);
t23 = cos(qJ(1));
t34 = t23 * t21;
t20 = sin(qJ(1));
t36 = t20 * t22;
t6 = t18 * t36 + t34;
t5 = t6 * t14;
t35 = t23 * t18;
t28 = (-(t13 * t21 - t14 * t18) * r_i_i_C(1) - (t13 * t18 + t14 * t21) * r_i_i_C(2)) * t19;
t24 = -t38 * t19 + t33 * t22;
t9 = t20 * t18 + t22 * t34;
t8 = -t20 * t21 + t22 * t35;
t7 = t21 * t36 - t35;
t2 = t8 * t13 + t9 * t14;
t1 = -t9 * t13 + t8 * t14;
t3 = [t23 * pkin(7) - t5 * r_i_i_C(2) - t39 * t20 - t26 * t7 + t27 * t6, t24 * t23, t25 * t9 - t26 * t8, t8, -t23 * t19, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t20 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t39 * t23 + t31 * t8 + t37 * t9, t24 * t20, -t5 * r_i_i_C(1) + t25 * t7 + t29 * t6, t6, -t20 * t19 (-t7 * t13 + t5) * r_i_i_C(1) + (-t6 * t13 - t7 * t14) * r_i_i_C(2); 0, t38 * t22 + t30 (-t37 * t18 + t31 * t21) * t19 - t28, t19 * t18, t22, t28;];
Ja_transl  = t3;
