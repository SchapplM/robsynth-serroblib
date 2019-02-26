% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:06
% EndTime: 2019-02-26 20:11:07
% DurationCPUTime: 0.16s
% Computational Cost: add. (216->32), mult. (342->55), div. (0->0), fcn. (426->12), ass. (0->32)
t52 = r_i_i_C(3) + qJ(5);
t27 = sin(pkin(12));
t30 = cos(pkin(12));
t51 = r_i_i_C(1) * t30 - r_i_i_C(2) * t27 + pkin(4);
t28 = sin(pkin(11));
t29 = sin(pkin(6));
t48 = t28 * t29;
t31 = cos(pkin(11));
t47 = t29 * t31;
t34 = sin(qJ(2));
t46 = t29 * t34;
t35 = cos(qJ(3));
t45 = t29 * t35;
t32 = cos(pkin(6));
t44 = t32 * t34;
t36 = cos(qJ(2));
t43 = t32 * t36;
t42 = t27 * r_i_i_C(1) + t30 * r_i_i_C(2) + pkin(8) + pkin(9);
t17 = t28 * t36 + t31 * t44;
t26 = qJ(3) + qJ(4);
t24 = sin(t26);
t25 = cos(t26);
t7 = t17 * t24 + t25 * t47;
t41 = t52 * (t17 * t25 - t24 * t47) - t51 * t7;
t19 = -t28 * t44 + t31 * t36;
t9 = t19 * t24 - t25 * t48;
t40 = t52 * (t19 * t25 + t24 * t48) - t51 * t9;
t14 = t24 * t46 - t32 * t25;
t39 = t52 * (t32 * t24 + t25 * t46) - t51 * t14;
t38 = t35 * pkin(3) + t52 * t24 + t51 * t25 + pkin(2);
t33 = sin(qJ(3));
t1 = [0, t42 * t19 + t38 * (-t28 * t43 - t31 * t34) (-t19 * t33 + t28 * t45) * pkin(3) + t40, t40, t9, 0; 0, t42 * t17 + t38 * (-t28 * t34 + t31 * t43) (-t17 * t33 - t31 * t45) * pkin(3) + t41, t41, t7, 0; 1 (t42 * t34 + t38 * t36) * t29 (t32 * t35 - t33 * t46) * pkin(3) + t39, t39, t14, 0;];
Ja_transl  = t1;
