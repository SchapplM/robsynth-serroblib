% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR8_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR8_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:51:32
% EndTime: 2019-02-26 22:51:32
% DurationCPUTime: 0.25s
% Computational Cost: add. (379->72), mult. (871->123), div. (0->0), fcn. (1119->14), ass. (0->52)
t43 = sin(qJ(4));
t75 = pkin(4) * t43 + pkin(10);
t47 = cos(qJ(4));
t36 = t47 * pkin(4) + pkin(3);
t39 = qJ(4) + qJ(5);
t37 = sin(t39);
t38 = cos(t39);
t54 = t38 * r_i_i_C(1) - t37 * r_i_i_C(2) + t36;
t74 = r_i_i_C(1) * t37 + r_i_i_C(2) * t38 + t75;
t45 = sin(qJ(2));
t46 = sin(qJ(1));
t49 = cos(qJ(2));
t50 = cos(qJ(1));
t60 = cos(pkin(6));
t56 = t50 * t60;
t28 = t46 * t45 - t49 * t56;
t29 = t45 * t56 + t46 * t49;
t44 = sin(qJ(3));
t48 = cos(qJ(3));
t40 = sin(pkin(7));
t41 = sin(pkin(6));
t67 = t41 * t50;
t58 = t40 * t67;
t42 = cos(pkin(7));
t66 = t42 * t44;
t12 = t28 * t66 - t29 * t48 + t44 * t58;
t21 = -t28 * t40 + t42 * t67;
t73 = (t12 * t37 - t21 * t38) * r_i_i_C(1) + (t12 * t38 + t21 * t37) * r_i_i_C(2);
t57 = t46 * t60;
t30 = -t50 * t45 - t49 * t57;
t31 = -t45 * t57 + t50 * t49;
t68 = t41 * t46;
t59 = t40 * t68;
t14 = t31 * t48 + (t30 * t42 + t59) * t44;
t23 = -t30 * t40 + t42 * t68;
t5 = -t14 * t37 + t23 * t38;
t6 = t14 * t38 + t23 * t37;
t72 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t55 = t60 * t40;
t62 = t45 * t48;
t63 = t44 * t49;
t20 = t44 * t55 + (t42 * t63 + t62) * t41;
t27 = -t41 * t49 * t40 + t60 * t42;
t71 = (-t20 * t37 + t27 * t38) * r_i_i_C(1) + (-t20 * t38 - t27 * t37) * r_i_i_C(2);
t69 = r_i_i_C(3) + pkin(12) + pkin(11);
t65 = t42 * t48;
t64 = t44 * t45;
t61 = t48 * t49;
t53 = t74 * t40;
t52 = -t29 * t44 + (-t28 * t42 - t58) * t48;
t13 = -t30 * t65 + t31 * t44 - t48 * t59;
t1 = [-t46 * pkin(1) - t29 * pkin(2) + pkin(9) * t67 + t54 * t12 + t74 * t21 + t69 * t52, t30 * pkin(2) + t54 * (t30 * t48 - t31 * t66) + t31 * t53 + t69 * (t30 * t44 + t31 * t65) -t54 * t13 + t69 * t14 (-t14 * t43 + t23 * t47) * pkin(4) + t72, t72, 0; t50 * pkin(1) + t31 * pkin(2) + pkin(9) * t68 + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t69 * t13 + t14 * t36 + t75 * t23, -t28 * pkin(2) + t54 * (-t28 * t48 - t29 * t66) + t29 * t53 + t69 * (-t28 * t44 + t29 * t65) -t12 * t69 + t54 * t52 (t12 * t43 - t21 * t47) * pkin(4) + t73, t73, 0; 0 (t69 * (t42 * t62 + t63) + t54 * (-t42 * t64 + t61) + t49 * pkin(2) + t45 * t53) * t41, t69 * t20 + t54 * (t48 * t55 + (t42 * t61 - t64) * t41) (-t20 * t43 + t27 * t47) * pkin(4) + t71, t71, 0;];
Ja_transl  = t1;
