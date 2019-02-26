% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:33:03
% EndTime: 2019-02-26 22:33:03
% DurationCPUTime: 0.20s
% Computational Cost: add. (280->52), mult. (380->72), div. (0->0), fcn. (449->10), ass. (0->43)
t31 = sin(qJ(4));
t35 = cos(qJ(4));
t62 = pkin(4) + pkin(5);
t67 = -qJ(5) * t31 - t62 * t35;
t29 = qJ(2) + qJ(3);
t27 = cos(t29);
t26 = sin(t29);
t30 = sin(qJ(6));
t34 = cos(qJ(6));
t45 = t30 * t31 + t34 * t35;
t64 = t45 * t26;
t46 = t30 * t35 - t31 * t34;
t9 = t46 * t26;
t66 = pkin(9) * t27 - r_i_i_C(1) * t64 + r_i_i_C(2) * t9;
t61 = -pkin(10) - r_i_i_C(3);
t65 = t27 * pkin(3) + (pkin(9) + t61) * t26;
t28 = cos(qJ(2)) * pkin(2);
t63 = pkin(1) + t28 + t65;
t36 = cos(qJ(1));
t55 = t36 * t35;
t33 = sin(qJ(1));
t58 = t33 * t31;
t12 = t27 * t58 + t55;
t11 = t12 * t34;
t57 = t33 * t35;
t56 = t36 * t31;
t52 = t66 * t33;
t51 = t66 * t36;
t49 = -t30 * r_i_i_C(1) - qJ(5);
t48 = t30 * r_i_i_C(2) - t62;
t47 = -t9 * r_i_i_C(1) - r_i_i_C(2) * t64;
t44 = t34 * r_i_i_C(2) - t49;
t43 = -t34 * r_i_i_C(1) + t48;
t40 = t65 + (r_i_i_C(1) * t45 - r_i_i_C(2) * t46 - t67) * t27;
t39 = t61 * t27 + (-pkin(3) + t67) * t26;
t38 = -sin(qJ(2)) * pkin(2) + t39;
t37 = -pkin(8) - pkin(7);
t15 = t27 * t55 + t58;
t14 = t27 * t56 - t57;
t13 = t27 * t57 - t56;
t2 = t14 * t30 + t15 * t34;
t1 = t14 * t34 - t15 * t30;
t3 = [-t11 * r_i_i_C(2) + t49 * t12 + t43 * t13 - t63 * t33 - t36 * t37, t36 * t38 + t51, t36 * t39 + t51, t43 * t14 + t44 * t15, t14, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t14 * qJ(5) + t62 * t15 - t33 * t37 + t63 * t36, t33 * t38 + t52, t33 * t39 + t52, -t11 * r_i_i_C(1) + t48 * t12 + t44 * t13, t12 (-t13 * t30 + t11) * r_i_i_C(1) + (-t12 * t30 - t13 * t34) * r_i_i_C(2); 0, t28 + t40, t40 (qJ(5) * t35 - t62 * t31) * t26 - t47, t26 * t31, t47;];
Ja_transl  = t3;
