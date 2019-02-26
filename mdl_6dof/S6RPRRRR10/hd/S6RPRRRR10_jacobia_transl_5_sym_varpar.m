% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR10_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR10_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:49
% EndTime: 2019-02-26 21:19:49
% DurationCPUTime: 0.18s
% Computational Cost: add. (290->57), mult. (669->95), div. (0->0), fcn. (869->14), ass. (0->48)
t38 = sin(qJ(4));
t66 = t38 * pkin(4) + pkin(9);
t33 = sin(pkin(7));
t36 = cos(pkin(7));
t37 = cos(pkin(6));
t32 = sin(pkin(13));
t43 = cos(qJ(1));
t54 = t43 * t32;
t35 = cos(pkin(13));
t40 = sin(qJ(1));
t55 = t40 * t35;
t47 = t37 * t55 + t54;
t34 = sin(pkin(6));
t56 = t40 * t34;
t65 = -t33 * t56 + t47 * t36;
t41 = cos(qJ(4));
t28 = t41 * pkin(4) + pkin(3);
t31 = qJ(4) + qJ(5);
t29 = sin(t31);
t30 = cos(t31);
t48 = t30 * r_i_i_C(1) - t29 * r_i_i_C(2) + t28;
t52 = t43 * t35;
t57 = t40 * t32;
t22 = -t37 * t52 + t57;
t23 = t37 * t54 + t55;
t39 = sin(qJ(3));
t42 = cos(qJ(3));
t53 = t43 * t34;
t50 = t33 * t53;
t58 = t36 * t39;
t12 = t22 * t58 - t23 * t42 + t39 * t50;
t17 = -t22 * t33 + t36 * t53;
t64 = (t12 * t29 - t17 * t30) * r_i_i_C(1) + (t12 * t30 + t17 * t29) * r_i_i_C(2);
t24 = -t37 * t57 + t52;
t14 = t24 * t42 - t65 * t39;
t19 = t47 * t33 + t36 * t56;
t5 = -t14 * t29 + t19 * t30;
t6 = t14 * t30 + t19 * t29;
t63 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t59 = t33 * t37;
t16 = t39 * t59 + (t32 * t42 + t35 * t58) * t34;
t21 = -t34 * t35 * t33 + t37 * t36;
t62 = (-t16 * t29 + t21 * t30) * r_i_i_C(1) + (-t16 * t30 - t21 * t29) * r_i_i_C(2);
t60 = r_i_i_C(3) + pkin(11) + pkin(10);
t51 = t34 * qJ(2);
t45 = -t23 * t39 + (-t22 * t36 - t50) * t42;
t13 = t24 * t39 + t65 * t42;
t1 = [-t40 * pkin(1) - t23 * pkin(2) + t43 * t51 + t60 * t45 + t48 * t12 + (t29 * r_i_i_C(1) + t30 * r_i_i_C(2) + t66) * t17, t56, -t48 * t13 + t60 * t14 (-t14 * t38 + t19 * t41) * pkin(4) + t63, t63, 0; t43 * pkin(1) + t24 * pkin(2) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t60 * t13 + t14 * t28 + t66 * t19 + t40 * t51, -t53, -t12 * t60 + t48 * t45 (t12 * t38 - t17 * t41) * pkin(4) + t64, t64, 0; 0, t37, t60 * t16 + t48 * (t42 * t59 + (t35 * t36 * t42 - t32 * t39) * t34) (-t16 * t38 + t21 * t41) * pkin(4) + t62, t62, 0;];
Ja_transl  = t1;
