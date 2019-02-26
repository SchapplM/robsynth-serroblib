% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP7
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:13
% EndTime: 2019-02-26 22:43:13
% DurationCPUTime: 0.20s
% Computational Cost: add. (375->57), mult. (592->91), div. (0->0), fcn. (739->12), ass. (0->44)
t64 = r_i_i_C(3) + qJ(6) + pkin(11);
t42 = cos(qJ(5));
t31 = t42 * pkin(5) + pkin(4);
t38 = sin(qJ(5));
t53 = t42 * r_i_i_C(1) - t38 * r_i_i_C(2) + t31;
t67 = pkin(5) + r_i_i_C(1);
t46 = -pkin(10) - pkin(9);
t49 = t42 * r_i_i_C(2) + t67 * t38 - t46;
t43 = cos(qJ(3));
t32 = t43 * pkin(3) + pkin(2);
t35 = qJ(3) + qJ(4);
t33 = sin(t35);
t34 = cos(t35);
t68 = t64 * t33 + t53 * t34 + t32;
t36 = sin(pkin(6));
t40 = sin(qJ(2));
t63 = t36 * t40;
t41 = sin(qJ(1));
t62 = t36 * t41;
t44 = cos(qJ(2));
t61 = t36 * t44;
t45 = cos(qJ(1));
t60 = t36 * t45;
t59 = cos(pkin(6));
t56 = t45 * t59;
t24 = t40 * t56 + t41 * t44;
t12 = t24 * t34 - t33 * t60;
t57 = t41 * t59;
t39 = sin(qJ(3));
t55 = t36 * (pkin(3) * t39 + pkin(8));
t26 = -t40 * t57 + t45 * t44;
t16 = t26 * t34 + t33 * t62;
t25 = t45 * t40 + t44 * t57;
t1 = -t16 * t38 + t25 * t42;
t11 = t24 * t33 + t34 * t60;
t51 = -t53 * t11 + t64 * t12;
t15 = t26 * t33 - t34 * t62;
t50 = -t53 * t15 + t64 * t16;
t21 = t33 * t63 - t59 * t34;
t22 = t59 * t33 + t34 * t63;
t48 = -t53 * t21 + t64 * t22;
t23 = t41 * t40 - t44 * t56;
t2 = t16 * t42 + t25 * t38;
t3 = [-t41 * pkin(1) - t64 * t11 - t53 * t12 - t49 * t23 - t24 * t32 + t45 * t55, -t25 * t68 + t49 * t26 (-t26 * t39 + t43 * t62) * pkin(3) + t50, t50, -t2 * r_i_i_C(2) + t67 * t1, t15; t45 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * t31 + t26 * t32 + t41 * t55 + (t38 * pkin(5) - t46) * t25 + t64 * t15, -t23 * t68 + t49 * t24 (-t24 * t39 - t43 * t60) * pkin(3) + t51, t51 (-t12 * t42 - t23 * t38) * r_i_i_C(2) + t67 * (-t12 * t38 + t23 * t42) t11; 0 (t49 * t40 + t68 * t44) * t36 (-t39 * t63 + t59 * t43) * pkin(3) + t48, t48 (-t22 * t42 + t38 * t61) * r_i_i_C(2) + t67 * (-t22 * t38 - t42 * t61) t21;];
Ja_transl  = t3;
