% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR15_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR15_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:25
% EndTime: 2019-02-26 22:24:25
% DurationCPUTime: 0.22s
% Computational Cost: add. (296->65), mult. (807->113), div. (0->0), fcn. (1042->12), ass. (0->47)
t65 = pkin(4) + pkin(10);
t37 = sin(qJ(2));
t38 = sin(qJ(1));
t41 = cos(qJ(2));
t42 = cos(qJ(1));
t52 = cos(pkin(6));
t47 = t42 * t52;
t24 = t37 * t47 + t38 * t41;
t36 = sin(qJ(3));
t40 = cos(qJ(3));
t23 = t38 * t37 - t41 * t47;
t34 = cos(pkin(7));
t32 = sin(pkin(7));
t33 = sin(pkin(6));
t59 = t33 * t42;
t49 = t32 * t59;
t44 = t23 * t34 + t49;
t64 = -t24 * t40 + t44 * t36;
t35 = sin(qJ(5));
t39 = cos(qJ(5));
t45 = t35 * r_i_i_C(1) + t39 * r_i_i_C(2) + qJ(4);
t63 = t39 * r_i_i_C(1) - t35 * r_i_i_C(2) + t65;
t62 = t24 * t36;
t60 = t33 * t38;
t58 = t34 * t36;
t57 = t34 * t40;
t56 = t36 * t37;
t55 = t36 * t41;
t54 = t37 * t40;
t53 = t40 * t41;
t51 = r_i_i_C(3) + pkin(11) + pkin(3);
t50 = t32 * t60;
t48 = t38 * t52;
t46 = t52 * t32;
t15 = -t23 * t32 + t34 * t59;
t25 = -t42 * t37 - t41 * t48;
t17 = -t25 * t32 + t34 * t60;
t43 = t63 * t32;
t26 = -t37 * t48 + t42 * t41;
t22 = -t33 * t41 * t32 + t52 * t34;
t13 = -t40 * t46 + (-t34 * t53 + t56) * t33;
t8 = t26 * t40 + (t25 * t34 + t50) * t36;
t7 = -t25 * t57 + t26 * t36 - t40 * t50;
t3 = t23 * t57 + t40 * t49 + t62;
t2 = t17 * t39 + t7 * t35;
t1 = -t17 * t35 + t7 * t39;
t4 = [-t24 * pkin(2) - t38 * pkin(1) + pkin(9) * t59 + t51 * t64 + t45 * (-t44 * t40 - t62) + t63 * t15, t25 * pkin(2) + t45 * (t25 * t36 + t26 * t57) + t26 * t43 + t51 * (t25 * t40 - t26 * t58) t45 * t8 - t51 * t7, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t42 * pkin(1) + t26 * pkin(2) + pkin(9) * t60 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(4) + t65 * t17 + t51 * t8, -t23 * pkin(2) + t45 * (-t23 * t36 + t24 * t57) + t24 * t43 + t51 * (-t23 * t40 - t24 * t58) -t51 * t3 - t45 * t64, t3 (t15 * t35 + t3 * t39) * r_i_i_C(1) + (t15 * t39 - t3 * t35) * r_i_i_C(2), 0; 0 (t45 * (t34 * t54 + t55) + t41 * pkin(2) + t37 * t43 + t51 * (-t34 * t56 + t53)) * t33, t45 * (t36 * t46 + (t34 * t55 + t54) * t33) - t51 * t13, t13 (t13 * t39 - t22 * t35) * r_i_i_C(1) + (-t13 * t35 - t22 * t39) * r_i_i_C(2), 0;];
Ja_transl  = t4;
