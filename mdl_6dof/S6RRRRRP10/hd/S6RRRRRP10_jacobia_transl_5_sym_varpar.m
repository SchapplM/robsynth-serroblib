% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP10_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP10_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:00
% EndTime: 2019-02-26 22:45:00
% DurationCPUTime: 0.20s
% Computational Cost: add. (261->53), mult. (508->89), div. (0->0), fcn. (640->12), ass. (0->40)
t30 = sin(qJ(3));
t34 = cos(qJ(3));
t33 = cos(qJ(4));
t24 = pkin(4) * t33 + pkin(3);
t27 = qJ(4) + qJ(5);
t25 = sin(t27);
t26 = cos(t27);
t40 = r_i_i_C(1) * t26 - r_i_i_C(2) * t25 + t24;
t50 = r_i_i_C(3) + pkin(11) + pkin(10);
t54 = t50 * t30 + t40 * t34 + pkin(2);
t31 = sin(qJ(2));
t32 = sin(qJ(1));
t35 = cos(qJ(2));
t45 = cos(pkin(6));
t49 = cos(qJ(1));
t41 = t45 * t49;
t18 = t31 * t41 + t32 * t35;
t28 = sin(pkin(6));
t43 = t28 * t49;
t10 = t18 * t34 - t30 * t43;
t17 = t32 * t31 - t35 * t41;
t53 = (-t10 * t25 + t17 * t26) * r_i_i_C(1) + (-t10 * t26 - t17 * t25) * r_i_i_C(2);
t42 = t32 * t45;
t20 = -t31 * t42 + t49 * t35;
t48 = t28 * t32;
t14 = t20 * t34 + t30 * t48;
t19 = t49 * t31 + t35 * t42;
t5 = -t14 * t25 + t19 * t26;
t6 = t14 * t26 + t19 * t25;
t52 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t47 = t28 * t34;
t16 = t45 * t30 + t31 * t47;
t46 = t28 * t35;
t51 = (-t16 * t25 - t26 * t46) * r_i_i_C(1) + (-t16 * t26 + t25 * t46) * r_i_i_C(2);
t29 = sin(qJ(4));
t44 = pkin(4) * t29 + pkin(9);
t39 = -t18 * t30 - t34 * t43;
t38 = r_i_i_C(1) * t25 + r_i_i_C(2) * t26 + t44;
t13 = t20 * t30 - t32 * t47;
t1 = [-t32 * pkin(1) - t18 * pkin(2) + pkin(8) * t43 - t40 * t10 - t38 * t17 + t50 * t39, -t19 * t54 + t38 * t20, -t40 * t13 + t50 * t14 (-t14 * t29 + t19 * t33) * pkin(4) + t52, t52, 0; t49 * pkin(1) + t20 * pkin(2) + pkin(8) * t48 + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t50 * t13 + t14 * t24 + t44 * t19, -t17 * t54 + t38 * t18, t50 * t10 + t40 * t39 (-t10 * t29 + t17 * t33) * pkin(4) + t53, t53, 0; 0 (t38 * t31 + t54 * t35) * t28, t50 * t16 + t40 * (-t28 * t31 * t30 + t45 * t34) (-t16 * t29 - t33 * t46) * pkin(4) + t51, t51, 0;];
Ja_transl  = t1;
