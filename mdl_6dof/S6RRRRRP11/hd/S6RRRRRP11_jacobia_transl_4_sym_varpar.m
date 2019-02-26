% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP11_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP11_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_jacobia_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:38
% EndTime: 2019-02-26 22:45:39
% DurationCPUTime: 0.25s
% Computational Cost: add. (242->73), mult. (663->129), div. (0->0), fcn. (854->12), ass. (0->49)
t35 = sin(qJ(2));
t36 = sin(qJ(1));
t39 = cos(qJ(2));
t40 = cos(qJ(1));
t48 = cos(pkin(6));
t44 = t40 * t48;
t22 = t36 * t35 - t39 * t44;
t30 = sin(pkin(7));
t32 = cos(pkin(7));
t31 = sin(pkin(6));
t55 = t31 * t40;
t15 = -t22 * t30 + t32 * t55;
t33 = sin(qJ(4));
t37 = cos(qJ(4));
t23 = t35 * t44 + t36 * t39;
t34 = sin(qJ(3));
t38 = cos(qJ(3));
t46 = t30 * t55;
t54 = t32 * t34;
t6 = t22 * t54 - t23 * t38 + t34 * t46;
t62 = t15 * t37 - t6 * t33;
t61 = t15 * t33 + t6 * t37;
t60 = r_i_i_C(3) + pkin(11);
t59 = t30 * pkin(10);
t58 = t30 * t33;
t57 = t30 * t37;
t56 = t31 * t36;
t53 = t32 * t38;
t52 = t34 * t35;
t51 = t34 * t39;
t50 = t35 * t38;
t49 = t38 * t39;
t47 = t30 * t56;
t45 = t36 * t48;
t43 = t48 * t30;
t42 = t37 * r_i_i_C(1) - t33 * r_i_i_C(2) + pkin(3);
t24 = -t40 * t35 - t39 * t45;
t17 = -t24 * t30 + t32 * t56;
t41 = -t23 * t34 + (-t22 * t32 - t46) * t38;
t25 = -t35 * t45 + t40 * t39;
t21 = -t31 * t39 * t30 + t48 * t32;
t14 = t34 * t43 + (t32 * t51 + t50) * t31;
t12 = t24 * t38 - t25 * t54;
t10 = -t22 * t38 - t23 * t54;
t8 = t25 * t38 + (t24 * t32 + t47) * t34;
t7 = -t24 * t53 + t25 * t34 - t38 * t47;
t2 = t17 * t33 + t8 * t37;
t1 = t17 * t37 - t8 * t33;
t3 = [-t36 * pkin(1) - t23 * pkin(2) + t6 * pkin(3) + pkin(9) * t55 + t15 * pkin(10) + t61 * r_i_i_C(1) + t62 * r_i_i_C(2) + t60 * t41 (t12 * t37 + t25 * t58) * r_i_i_C(1) + (-t12 * t33 + t25 * t57) * r_i_i_C(2) + t12 * pkin(3) + t24 * pkin(2) + t25 * t59 + t60 * (t24 * t34 + t25 * t53) -t42 * t7 + t60 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t40 * pkin(1) + t25 * pkin(2) + t8 * pkin(3) + pkin(9) * t56 + t17 * pkin(10) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t60 * t7 (t10 * t37 + t23 * t58) * r_i_i_C(1) + (-t10 * t33 + t23 * t57) * r_i_i_C(2) + t10 * pkin(3) - t22 * pkin(2) + t23 * t59 + t60 * (-t22 * t34 + t23 * t53) t42 * t41 - t6 * t60, -t62 * r_i_i_C(1) + t61 * r_i_i_C(2), 0, 0; 0 (t42 * (-t32 * t52 + t49) + t39 * pkin(2) + (t33 * r_i_i_C(1) + t37 * r_i_i_C(2) + pkin(10)) * t35 * t30 + t60 * (t32 * t50 + t51)) * t31, t60 * t14 + t42 * (t38 * t43 + (t32 * t49 - t52) * t31) (-t14 * t33 + t21 * t37) * r_i_i_C(1) + (-t14 * t37 - t21 * t33) * r_i_i_C(2), 0, 0;];
Ja_transl  = t3;
