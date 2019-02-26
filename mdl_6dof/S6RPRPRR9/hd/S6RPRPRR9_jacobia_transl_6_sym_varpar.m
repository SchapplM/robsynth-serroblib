% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:23
% EndTime: 2019-02-26 20:53:24
% DurationCPUTime: 0.38s
% Computational Cost: add. (658->83), mult. (1787->137), div. (0->0), fcn. (2376->16), ass. (0->60)
t42 = sin(pkin(13));
t46 = cos(pkin(13));
t52 = sin(qJ(3));
t56 = cos(qJ(3));
t37 = t52 * t42 - t56 * t46;
t44 = sin(pkin(7));
t26 = t37 * t44;
t48 = cos(pkin(7));
t28 = t37 * t48;
t49 = cos(pkin(6));
t47 = cos(pkin(12));
t57 = cos(qJ(1));
t66 = t57 * t47;
t43 = sin(pkin(12));
t53 = sin(qJ(1));
t71 = t53 * t43;
t33 = -t49 * t66 + t71;
t68 = t57 * t43;
t69 = t53 * t47;
t34 = t49 * t68 + t69;
t62 = t56 * t42 + t52 * t46;
t45 = sin(pkin(6));
t67 = t57 * t45;
t14 = -t26 * t67 - t33 * t28 + t34 * t62;
t50 = sin(qJ(6));
t54 = cos(qJ(6));
t27 = t62 * t44;
t29 = t62 * t48;
t15 = t27 * t67 + t33 * t29 + t34 * t37;
t23 = -t33 * t44 + t48 * t67;
t51 = sin(qJ(5));
t55 = cos(qJ(5));
t6 = t15 * t55 + t23 * t51;
t79 = t14 * t54 + t6 * t50;
t78 = -t14 * t50 + t6 * t54;
t75 = t15 * t51 - t23 * t55;
t21 = t49 * t27 + (t29 * t47 - t37 * t43) * t45;
t74 = pkin(11) + r_i_i_C(3);
t73 = pkin(3) * t52;
t70 = t53 * t45;
t65 = pkin(9) + qJ(4);
t64 = t45 * (t44 * t73 + t65 * t48 + qJ(2));
t61 = t54 * r_i_i_C(1) - t50 * r_i_i_C(2) + pkin(5);
t60 = -t50 * r_i_i_C(1) - t54 * r_i_i_C(2) - pkin(10);
t35 = -t49 * t69 - t68;
t59 = -t35 * t44 + t48 * t70;
t36 = -t49 * t71 + t66;
t18 = t27 * t70 + t35 * t29 - t36 * t37;
t58 = t74 * t51 + t61 * t55 + pkin(4);
t41 = t56 * pkin(3) + pkin(2);
t32 = -t45 * t47 * t44 + t49 * t48;
t31 = -t65 * t44 + t48 * t73;
t20 = -t49 * t26 + (-t28 * t47 - t43 * t62) * t45;
t17 = -t26 * t70 - t35 * t28 - t36 * t62;
t10 = t21 * t55 + t32 * t51;
t8 = t18 * t55 + t59 * t51;
t7 = t18 * t51 - t59 * t55;
t2 = -t17 * t50 + t8 * t54;
t1 = -t17 * t54 - t8 * t50;
t3 = [-t53 * pkin(1) + t15 * pkin(4) + t6 * pkin(5) - t14 * pkin(10) + t78 * r_i_i_C(1) - t79 * r_i_i_C(2) + t33 * t31 - t34 * t41 + t57 * t64 + t74 * t75, t70, -t60 * t18 + (-t36 * t52 + (t35 * t48 + t44 * t70) * t56) * pkin(3) + t58 * t17, t59, -t61 * t7 + t74 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t57 * pkin(1) + t18 * pkin(4) + t8 * pkin(5) - t17 * pkin(10) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t35 * t31 + t36 * t41 + t53 * t64 + t74 * t7, -t67, t60 * t15 + (-t34 * t52 + (-t33 * t48 - t44 * t67) * t56) * pkin(3) - t58 * t14, -t23, -t6 * t74 + t61 * t75, t79 * r_i_i_C(1) + t78 * r_i_i_C(2); 0, t49, -t60 * t21 + (t49 * t44 * t56 + (t47 * t48 * t56 - t43 * t52) * t45) * pkin(3) + t58 * t20, t32, t74 * t10 + t61 * (-t21 * t51 + t32 * t55) (-t10 * t50 - t20 * t54) * r_i_i_C(1) + (-t10 * t54 + t20 * t50) * r_i_i_C(2);];
Ja_transl  = t3;
