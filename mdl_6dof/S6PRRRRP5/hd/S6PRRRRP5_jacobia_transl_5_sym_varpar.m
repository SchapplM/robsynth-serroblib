% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRP5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRP5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:17:28
% EndTime: 2019-02-26 20:17:29
% DurationCPUTime: 0.29s
% Computational Cost: add. (403->84), mult. (1135->157), div. (0->0), fcn. (1481->14), ass. (0->59)
t36 = sin(pkin(7));
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t38 = cos(pkin(12));
t64 = cos(pkin(6));
t58 = t38 * t64;
t62 = sin(pkin(12));
t49 = t62 * t42 - t45 * t58;
t63 = cos(pkin(7));
t37 = sin(pkin(6));
t67 = t37 * t38;
t74 = t36 * t67 + t49 * t63;
t53 = t64 * t62;
t50 = t38 * t42 + t45 * t53;
t59 = t37 * t62;
t73 = -t36 * t59 + t50 * t63;
t72 = r_i_i_C(3) + pkin(11);
t71 = pkin(9) * t36;
t70 = cos(qJ(3));
t40 = sin(qJ(4));
t69 = t36 * t40;
t44 = cos(qJ(4));
t68 = t36 * t44;
t66 = t37 * t42;
t65 = t37 * t45;
t61 = t36 * t66;
t41 = sin(qJ(3));
t57 = t41 * t63;
t56 = t64 * t36;
t54 = t63 * t70;
t39 = sin(qJ(5));
t43 = cos(qJ(5));
t52 = t43 * r_i_i_C(1) - t39 * r_i_i_C(2) + pkin(4);
t51 = t39 * r_i_i_C(1) + t43 * r_i_i_C(2) + pkin(10);
t48 = -t72 * t40 - t52 * t44 - pkin(3);
t31 = t38 * t45 - t42 * t53;
t30 = t42 * t58 + t62 * t45;
t29 = -t36 * t65 + t64 * t63;
t28 = (-t42 * t57 + t70 * t45) * t37;
t27 = (t41 * t45 + t42 * t54) * t37;
t24 = t50 * t36 + t63 * t59;
t23 = t49 * t36 - t63 * t67;
t22 = t41 * t56 + (t70 * t42 + t45 * t57) * t37;
t21 = t41 * t66 - t54 * t65 - t70 * t56;
t20 = t28 * t44 + t40 * t61;
t18 = -t31 * t57 - t50 * t70;
t17 = t31 * t54 - t50 * t41;
t16 = -t30 * t57 - t49 * t70;
t15 = t30 * t54 - t49 * t41;
t14 = t22 * t44 + t29 * t40;
t12 = t31 * t70 - t73 * t41;
t11 = t31 * t41 + t73 * t70;
t10 = t30 * t70 - t74 * t41;
t9 = t30 * t41 + t74 * t70;
t8 = t18 * t44 + t31 * t69;
t6 = t16 * t44 + t30 * t69;
t4 = t12 * t44 + t24 * t40;
t2 = t10 * t44 + t23 * t40;
t1 = [0 (t17 * t39 + t8 * t43) * r_i_i_C(1) + (t17 * t43 - t8 * t39) * r_i_i_C(2) + t8 * pkin(4) + t18 * pkin(3) + t17 * pkin(10) - t50 * pkin(2) + t31 * t71 + t72 * (t18 * t40 - t31 * t68) t48 * t11 + t51 * t12, t72 * t4 + t52 * (-t12 * t40 + t24 * t44) (t11 * t43 - t4 * t39) * r_i_i_C(1) + (-t11 * t39 - t4 * t43) * r_i_i_C(2), 0; 0 (t15 * t39 + t6 * t43) * r_i_i_C(1) + (t15 * t43 - t6 * t39) * r_i_i_C(2) + t6 * pkin(4) + t16 * pkin(3) + t15 * pkin(10) - t49 * pkin(2) + t30 * t71 + t72 * (t16 * t40 - t30 * t68) t51 * t10 + t48 * t9, t72 * t2 + t52 * (-t10 * t40 + t23 * t44) (-t2 * t39 + t9 * t43) * r_i_i_C(1) + (-t2 * t43 - t9 * t39) * r_i_i_C(2), 0; 1 (t20 * t43 + t27 * t39) * r_i_i_C(1) + (-t20 * t39 + t27 * t43) * r_i_i_C(2) + t20 * pkin(4) + t28 * pkin(3) + t27 * pkin(10) + (t45 * pkin(2) + t42 * t71) * t37 + t72 * (t28 * t40 - t44 * t61) t48 * t21 + t51 * t22, t72 * t14 + t52 * (-t22 * t40 + t29 * t44) (-t14 * t39 + t21 * t43) * r_i_i_C(1) + (-t14 * t43 - t21 * t39) * r_i_i_C(2), 0;];
Ja_transl  = t1;
