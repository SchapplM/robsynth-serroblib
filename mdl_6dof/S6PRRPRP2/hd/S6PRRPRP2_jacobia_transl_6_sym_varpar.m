% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRP2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:47
% EndTime: 2019-02-26 20:01:47
% DurationCPUTime: 0.21s
% Computational Cost: add. (284->54), mult. (484->94), div. (0->0), fcn. (617->12), ass. (0->40)
t29 = qJ(3) + pkin(11);
t27 = sin(t29);
t28 = cos(t29);
t39 = cos(qJ(3));
t55 = pkin(9) + r_i_i_C(2);
t57 = t39 * pkin(3) + pkin(4) * t28 + t55 * t27 + pkin(2);
t56 = pkin(5) + r_i_i_C(1);
t35 = sin(qJ(5));
t54 = t28 * t35;
t38 = cos(qJ(5));
t53 = t28 * t38;
t30 = sin(pkin(10));
t31 = sin(pkin(6));
t52 = t30 * t31;
t32 = cos(pkin(10));
t51 = t31 * t32;
t37 = sin(qJ(2));
t50 = t31 * t37;
t49 = t31 * t39;
t33 = cos(pkin(6));
t48 = t33 * t37;
t40 = cos(qJ(2));
t47 = t33 * t40;
t46 = t35 * t40;
t45 = t38 * t40;
t44 = r_i_i_C(3) + qJ(6);
t41 = t44 * t35 + t56 * t38 + pkin(4);
t36 = sin(qJ(3));
t34 = -qJ(4) - pkin(8);
t24 = -t30 * t48 + t32 * t40;
t23 = t30 * t47 + t32 * t37;
t22 = t30 * t40 + t32 * t48;
t21 = t30 * t37 - t32 * t47;
t18 = t33 * t27 + t28 * t50;
t13 = t18 * t35 + t31 * t45;
t12 = t24 * t28 + t27 * t52;
t10 = t22 * t28 - t27 * t51;
t3 = t12 * t35 - t23 * t38;
t1 = t10 * t35 - t21 * t38;
t2 = [0, -t24 * t34 + t56 * (-t23 * t53 + t24 * t35) + t44 * (-t23 * t54 - t24 * t38) - t57 * t23, t55 * t12 + (-t24 * t36 + t30 * t49) * pkin(3) + t41 * (-t24 * t27 + t28 * t52) t23, t44 * (t12 * t38 + t23 * t35) - t56 * t3, t3; 0, -t22 * t34 + t56 * (-t21 * t53 + t22 * t35) + t44 * (-t21 * t54 - t22 * t38) - t57 * t21, t55 * t10 + (-t22 * t36 - t32 * t49) * pkin(3) + t41 * (-t22 * t27 - t28 * t51) t21, t44 * (t10 * t38 + t21 * t35) - t56 * t1, t1; 1 (t56 * (t28 * t45 + t35 * t37) + t44 * (t28 * t46 - t37 * t38) - t34 * t37 + t57 * t40) * t31, t55 * t18 + (t33 * t39 - t36 * t50) * pkin(3) + t41 * (-t27 * t50 + t33 * t28) -t31 * t40, t44 * (t18 * t38 - t31 * t46) - t56 * t13, t13;];
Ja_transl  = t2;
