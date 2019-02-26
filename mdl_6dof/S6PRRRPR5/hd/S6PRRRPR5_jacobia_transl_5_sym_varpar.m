% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:55
% EndTime: 2019-02-26 20:12:56
% DurationCPUTime: 0.19s
% Computational Cost: add. (246->60), mult. (582->109), div. (0->0), fcn. (747->14), ass. (0->46)
t59 = r_i_i_C(3) + qJ(5) + pkin(10);
t30 = sin(pkin(7));
t31 = sin(pkin(6));
t58 = t30 * t31;
t33 = cos(pkin(7));
t57 = t31 * t33;
t36 = sin(qJ(3));
t56 = t33 * t36;
t39 = cos(qJ(3));
t55 = t33 * t39;
t37 = sin(qJ(2));
t54 = t36 * t37;
t40 = cos(qJ(2));
t53 = t36 * t40;
t52 = t37 * t39;
t51 = t39 * t40;
t50 = cos(pkin(6));
t49 = sin(pkin(12));
t32 = cos(pkin(12));
t48 = t32 * t58;
t47 = t31 * t49;
t46 = t32 * t50;
t45 = t50 * t30;
t44 = t30 * t47;
t43 = t50 * t49;
t29 = qJ(4) + pkin(13);
t27 = sin(t29);
t28 = cos(t29);
t38 = cos(qJ(4));
t42 = t38 * pkin(4) + t28 * r_i_i_C(1) - t27 * r_i_i_C(2) + pkin(3);
t35 = sin(qJ(4));
t41 = (t35 * pkin(4) + t27 * r_i_i_C(1) + t28 * r_i_i_C(2) + pkin(9)) * t30;
t21 = t32 * t40 - t37 * t43;
t20 = -t32 * t37 - t40 * t43;
t19 = t37 * t46 + t49 * t40;
t18 = -t49 * t37 + t40 * t46;
t17 = t50 * t33 - t40 * t58;
t12 = -t20 * t30 + t33 * t47;
t11 = -t18 * t30 - t32 * t57;
t10 = t36 * t45 + (t33 * t53 + t52) * t31;
t9 = t31 * t54 - t39 * t45 - t51 * t57;
t4 = t21 * t39 + (t20 * t33 + t44) * t36;
t3 = -t20 * t55 + t21 * t36 - t39 * t44;
t2 = t19 * t39 + (t18 * t33 - t48) * t36;
t1 = -t18 * t55 + t19 * t36 + t39 * t48;
t5 = [0, t20 * pkin(2) + t59 * (t20 * t36 + t21 * t55) + t42 * (t20 * t39 - t21 * t56) + t21 * t41, -t42 * t3 + t59 * t4 (t12 * t28 - t4 * t27) * r_i_i_C(1) + (-t12 * t27 - t4 * t28) * r_i_i_C(2) + (t12 * t38 - t4 * t35) * pkin(4), t3, 0; 0, t18 * pkin(2) + t59 * (t18 * t36 + t19 * t55) + t42 * (t18 * t39 - t19 * t56) + t19 * t41, -t42 * t1 + t59 * t2 (t11 * t28 - t2 * t27) * r_i_i_C(1) + (-t11 * t27 - t2 * t28) * r_i_i_C(2) + (t11 * t38 - t2 * t35) * pkin(4), t1, 0; 1 (t59 * (t33 * t52 + t53) + t42 * (-t33 * t54 + t51) + t40 * pkin(2) + t37 * t41) * t31, t59 * t10 - t42 * t9 (-t10 * t27 + t17 * t28) * r_i_i_C(1) + (-t10 * t28 - t17 * t27) * r_i_i_C(2) + (-t10 * t35 + t17 * t38) * pkin(4), t9, 0;];
Ja_transl  = t5;
