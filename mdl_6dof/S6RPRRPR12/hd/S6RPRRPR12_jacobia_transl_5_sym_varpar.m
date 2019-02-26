% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR12_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR12_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:15
% EndTime: 2019-02-26 21:07:16
% DurationCPUTime: 0.21s
% Computational Cost: add. (294->51), mult. (800->84), div. (0->0), fcn. (1050->12), ass. (0->45)
t33 = cos(pkin(6));
t31 = cos(pkin(12));
t39 = cos(qJ(1));
t49 = t39 * t31;
t28 = sin(pkin(12));
t36 = sin(qJ(1));
t54 = t36 * t28;
t22 = -t33 * t49 + t54;
t51 = t39 * t28;
t52 = t36 * t31;
t23 = t33 * t51 + t52;
t35 = sin(qJ(3));
t38 = cos(qJ(3));
t29 = sin(pkin(7));
t30 = sin(pkin(6));
t50 = t39 * t30;
t46 = t29 * t50;
t32 = cos(pkin(7));
t55 = t32 * t35;
t12 = t22 * t55 - t23 * t38 + t35 * t46;
t18 = -t22 * t29 + t32 * t50;
t34 = sin(qJ(4));
t37 = cos(qJ(4));
t61 = t12 * t37 + t18 * t34;
t60 = t12 * t34 - t18 * t37;
t44 = t33 * t52 + t51;
t53 = t36 * t30;
t59 = -t29 * t53 + t44 * t32;
t48 = r_i_i_C(3) + qJ(5);
t57 = -r_i_i_C(2) + pkin(4);
t42 = t48 * t34 + t57 * t37 + pkin(3);
t58 = r_i_i_C(1) + pkin(10);
t56 = t29 * t33;
t47 = t30 * qJ(2);
t41 = -t23 * t35 + (-t22 * t32 - t46) * t38;
t40 = t44 * t29 + t32 * t53;
t24 = -t33 * t54 + t49;
t21 = -t30 * t31 * t29 + t33 * t32;
t17 = t35 * t56 + (t28 * t38 + t31 * t55) * t30;
t14 = t24 * t38 - t35 * t59;
t13 = t24 * t35 + t38 * t59;
t7 = t17 * t34 - t21 * t37;
t6 = t14 * t37 + t40 * t34;
t5 = t14 * t34 - t40 * t37;
t1 = [-t36 * pkin(1) - t23 * pkin(2) + t12 * pkin(3) + t18 * pkin(9) + t39 * t47 + t58 * t41 + t48 * t60 + t57 * t61, t53, -t13 * t42 + t58 * t14, t48 * t6 - t57 * t5, t5, 0; t39 * pkin(1) + t24 * pkin(2) + t14 * pkin(3) + t40 * pkin(9) + t58 * t13 + t36 * t47 + t48 * t5 + t57 * t6, -t50, -t12 * t58 + t42 * t41, -t48 * t61 + t57 * t60, -t60, 0; 0, t33, t58 * t17 + t42 * (t38 * t56 + (t31 * t32 * t38 - t28 * t35) * t30) t48 * (t17 * t37 + t21 * t34) - t57 * t7, t7, 0;];
Ja_transl  = t1;
