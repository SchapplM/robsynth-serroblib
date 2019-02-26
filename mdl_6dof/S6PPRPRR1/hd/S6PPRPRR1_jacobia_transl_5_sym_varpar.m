% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRPRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRPRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:47
% EndTime: 2019-02-26 19:39:48
% DurationCPUTime: 0.20s
% Computational Cost: add. (183->46), mult. (511->90), div. (0->0), fcn. (675->14), ass. (0->38)
t27 = sin(pkin(7));
t24 = sin(pkin(13));
t29 = cos(pkin(13));
t35 = sin(qJ(3));
t37 = cos(qJ(3));
t40 = t37 * t24 + t35 * t29;
t13 = t40 * t27;
t32 = cos(pkin(7));
t15 = t40 * t32;
t21 = t35 * t24 - t37 * t29;
t25 = sin(pkin(12));
t28 = sin(pkin(6));
t30 = cos(pkin(12));
t33 = cos(pkin(6));
t9 = t33 * t13 + (t15 * t30 - t21 * t25) * t28;
t48 = -pkin(9) - r_i_i_C(3);
t26 = sin(pkin(11));
t47 = t26 * t28;
t46 = t26 * t33;
t45 = t27 * t28;
t31 = cos(pkin(11));
t44 = t31 * t28;
t43 = t31 * t33;
t34 = sin(qJ(5));
t36 = cos(qJ(5));
t39 = t36 * r_i_i_C(1) - t34 * r_i_i_C(2) + pkin(4);
t17 = -t26 * t25 + t30 * t43;
t18 = t25 * t43 + t26 * t30;
t38 = t13 * t44 - t17 * t15 + t18 * t21;
t19 = -t31 * t25 - t30 * t46;
t20 = -t25 * t46 + t31 * t30;
t6 = t13 * t47 + t19 * t15 - t20 * t21;
t16 = -t30 * t45 + t33 * t32;
t14 = t21 * t32;
t12 = t21 * t27;
t11 = -t19 * t27 + t32 * t47;
t10 = -t17 * t27 - t32 * t44;
t1 = [0, t47, -t48 * t6 + t39 * (-t12 * t47 - t19 * t14 - t20 * t40) + (-t20 * t35 + (t19 * t32 + t26 * t45) * t37) * pkin(3), t11 (t11 * t36 - t6 * t34) * r_i_i_C(1) + (-t11 * t34 - t6 * t36) * r_i_i_C(2), 0; 0, -t44, t48 * t38 + t39 * (t12 * t44 - t17 * t14 - t18 * t40) + (-t18 * t35 + (t17 * t32 - t27 * t44) * t37) * pkin(3), t10 (t10 * t36 + t34 * t38) * r_i_i_C(1) + (-t10 * t34 + t36 * t38) * r_i_i_C(2), 0; 1, t33, -t48 * t9 + t39 * (-t33 * t12 + (-t14 * t30 - t25 * t40) * t28) + (t27 * t33 * t37 + (t30 * t32 * t37 - t25 * t35) * t28) * pkin(3), t16 (t16 * t36 - t9 * t34) * r_i_i_C(1) + (-t16 * t34 - t9 * t36) * r_i_i_C(2), 0;];
Ja_transl  = t1;
