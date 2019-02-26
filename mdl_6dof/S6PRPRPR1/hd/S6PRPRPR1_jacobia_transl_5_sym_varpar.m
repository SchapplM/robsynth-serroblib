% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:46:27
% EndTime: 2019-02-26 19:46:27
% DurationCPUTime: 0.13s
% Computational Cost: add. (120->37), mult. (254->63), div. (0->0), fcn. (329->12), ass. (0->29)
t19 = sin(pkin(11));
t22 = cos(pkin(11));
t27 = sin(qJ(2));
t29 = cos(qJ(2));
t11 = t27 * t19 - t29 * t22;
t38 = r_i_i_C(3) + qJ(5) + pkin(8);
t20 = sin(pkin(10));
t21 = sin(pkin(6));
t37 = t20 * t21;
t23 = cos(pkin(10));
t36 = t23 * t21;
t24 = cos(pkin(6));
t35 = t24 * t29;
t31 = t29 * t19 + t27 * t22;
t10 = t31 * t24;
t3 = t23 * t10 - t20 * t11;
t32 = t20 * t10 + t23 * t11;
t18 = qJ(4) + pkin(12);
t16 = sin(t18);
t17 = cos(t18);
t28 = cos(qJ(4));
t30 = t28 * pkin(4) + t17 * r_i_i_C(1) - t16 * r_i_i_C(2) + pkin(3);
t26 = sin(qJ(4));
t9 = t11 * t24;
t8 = t31 * t21;
t7 = t11 * t21;
t5 = t20 * t9 - t23 * t31;
t2 = -t20 * t31 - t23 * t9;
t1 = [0, -t38 * t32 + (-t20 * t35 - t23 * t27) * pkin(2) + t30 * t5, t37 (t16 * t32 + t17 * t37) * r_i_i_C(1) + (-t16 * t37 + t17 * t32) * r_i_i_C(2) + (t26 * t32 + t28 * t37) * pkin(4), -t5, 0; 0, t38 * t3 + (-t20 * t27 + t23 * t35) * pkin(2) + t30 * t2, -t36 (-t3 * t16 - t17 * t36) * r_i_i_C(1) + (t16 * t36 - t3 * t17) * r_i_i_C(2) + (-t3 * t26 - t28 * t36) * pkin(4), -t2, 0; 1, t21 * t29 * pkin(2) - t30 * t7 + t38 * t8, t24 (-t8 * t16 + t24 * t17) * r_i_i_C(1) + (-t24 * t16 - t8 * t17) * r_i_i_C(2) + (t24 * t28 - t8 * t26) * pkin(4), t7, 0;];
Ja_transl  = t1;
