% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR6_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR6_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:04
% EndTime: 2019-02-26 20:07:05
% DurationCPUTime: 0.16s
% Computational Cost: add. (225->54), mult. (524->100), div. (0->0), fcn. (674->14), ass. (0->44)
t58 = r_i_i_C(3) + pkin(10) + qJ(4);
t31 = sin(pkin(7));
t32 = sin(pkin(6));
t57 = t31 * t32;
t34 = cos(pkin(7));
t56 = t32 * t34;
t36 = sin(qJ(3));
t55 = t34 * t36;
t38 = cos(qJ(3));
t54 = t34 * t38;
t37 = sin(qJ(2));
t53 = t36 * t37;
t39 = cos(qJ(2));
t52 = t36 * t39;
t51 = t37 * t38;
t50 = t38 * t39;
t49 = cos(pkin(6));
t48 = sin(pkin(12));
t33 = cos(pkin(12));
t47 = t33 * t57;
t46 = t32 * t48;
t45 = t33 * t49;
t44 = t49 * t31;
t43 = t31 * t46;
t42 = t49 * t48;
t29 = pkin(13) + qJ(5);
t27 = sin(t29);
t28 = cos(t29);
t41 = t28 * r_i_i_C(1) - t27 * r_i_i_C(2) + cos(pkin(13)) * pkin(4) + pkin(3);
t40 = (pkin(4) * sin(pkin(13)) + t27 * r_i_i_C(1) + t28 * r_i_i_C(2) + pkin(9)) * t31;
t21 = t33 * t39 - t37 * t42;
t20 = -t33 * t37 - t39 * t42;
t19 = t37 * t45 + t48 * t39;
t18 = -t48 * t37 + t39 * t45;
t17 = t49 * t34 - t39 * t57;
t12 = -t20 * t31 + t34 * t46;
t11 = -t18 * t31 - t33 * t56;
t10 = t36 * t44 + (t34 * t52 + t51) * t32;
t9 = t32 * t53 - t38 * t44 - t50 * t56;
t4 = t21 * t38 + (t20 * t34 + t43) * t36;
t3 = -t20 * t54 + t21 * t36 - t38 * t43;
t2 = t19 * t38 + (t18 * t34 - t47) * t36;
t1 = -t18 * t54 + t19 * t36 + t38 * t47;
t5 = [0, t20 * pkin(2) + t58 * (t20 * t36 + t21 * t54) + t41 * (t20 * t38 - t21 * t55) + t21 * t40, -t41 * t3 + t58 * t4, t3 (t12 * t28 - t4 * t27) * r_i_i_C(1) + (-t12 * t27 - t4 * t28) * r_i_i_C(2), 0; 0, t18 * pkin(2) + t58 * (t18 * t36 + t19 * t54) + t41 * (t18 * t38 - t19 * t55) + t19 * t40, -t41 * t1 + t58 * t2, t1 (t11 * t28 - t2 * t27) * r_i_i_C(1) + (-t11 * t27 - t2 * t28) * r_i_i_C(2), 0; 1 (t58 * (t34 * t51 + t52) + t41 * (-t34 * t53 + t50) + t39 * pkin(2) + t37 * t40) * t32, t58 * t10 - t41 * t9, t9 (-t10 * t27 + t17 * t28) * r_i_i_C(1) + (-t10 * t28 - t17 * t27) * r_i_i_C(2), 0;];
Ja_transl  = t5;
