% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRP4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:00
% EndTime: 2019-02-26 19:52:01
% DurationCPUTime: 0.14s
% Computational Cost: add. (168->36), mult. (279->64), div. (0->0), fcn. (353->11), ass. (0->29)
t15 = pkin(11) + qJ(4);
t13 = sin(t15);
t14 = cos(t15);
t19 = sin(qJ(5));
t21 = cos(qJ(5));
t25 = t21 * r_i_i_C(1) - t19 * r_i_i_C(2) + pkin(4);
t34 = pkin(9) + r_i_i_C(3);
t35 = t34 * t13 + t25 * t14 + cos(pkin(11)) * pkin(3) + pkin(2);
t16 = sin(pkin(10));
t17 = sin(pkin(6));
t33 = t16 * t17;
t20 = sin(qJ(2));
t32 = t17 * t20;
t22 = cos(qJ(2));
t31 = t17 * t22;
t30 = cos(pkin(6));
t29 = cos(pkin(10));
t28 = t16 * t30;
t27 = t17 * t29;
t26 = t30 * t29;
t24 = t19 * r_i_i_C(1) + t21 * r_i_i_C(2) + pkin(8) + qJ(3);
t10 = -t20 * t28 + t29 * t22;
t9 = t29 * t20 + t22 * t28;
t8 = t16 * t22 + t20 * t26;
t7 = t16 * t20 - t22 * t26;
t6 = t30 * t13 + t14 * t32;
t4 = t10 * t14 + t13 * t33;
t2 = -t13 * t27 + t8 * t14;
t1 = [0, t24 * t10 - t35 * t9, t9, t34 * t4 + t25 * (-t10 * t13 + t14 * t33) (-t4 * t19 + t9 * t21) * r_i_i_C(1) + (-t9 * t19 - t4 * t21) * r_i_i_C(2), 0; 0, t24 * t8 - t35 * t7, t7, t34 * t2 + t25 * (-t8 * t13 - t14 * t27) (-t2 * t19 + t7 * t21) * r_i_i_C(1) + (-t7 * t19 - t2 * t21) * r_i_i_C(2), 0; 1 (t24 * t20 + t35 * t22) * t17, -t31, t34 * t6 + t25 * (-t13 * t32 + t30 * t14) (-t6 * t19 - t21 * t31) * r_i_i_C(1) + (t19 * t31 - t6 * t21) * r_i_i_C(2), 0;];
Ja_transl  = t1;
