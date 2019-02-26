% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:47
% EndTime: 2019-02-26 19:49:47
% DurationCPUTime: 0.12s
% Computational Cost: add. (133->36), mult. (342->63), div. (0->0), fcn. (435->10), ass. (0->28)
t15 = sin(pkin(6));
t19 = sin(qJ(4));
t33 = t15 * t19;
t20 = sin(qJ(2));
t32 = t15 * t20;
t22 = cos(qJ(4));
t31 = t15 * t22;
t23 = cos(qJ(2));
t30 = t15 * t23;
t17 = cos(pkin(6));
t29 = t17 * t20;
t28 = t17 * t23;
t27 = pkin(4) + pkin(9) + r_i_i_C(3);
t18 = sin(qJ(6));
t21 = cos(qJ(6));
t26 = t18 * r_i_i_C(1) + t21 * r_i_i_C(2) + qJ(5);
t25 = t21 * r_i_i_C(1) - t18 * r_i_i_C(2) + pkin(2) + pkin(5) + pkin(8);
t24 = t27 * t19 - t26 * t22 + qJ(3);
t16 = cos(pkin(10));
t14 = sin(pkin(10));
t11 = t17 * t19 + t22 * t30;
t10 = -t14 * t29 + t16 * t23;
t9 = t14 * t28 + t16 * t20;
t8 = t14 * t23 + t16 * t29;
t7 = t14 * t20 - t16 * t28;
t3 = t16 * t33 + t7 * t22;
t1 = t14 * t33 - t9 * t22;
t2 = [0, t24 * t10 - t25 * t9, t9, t26 * (t14 * t31 + t9 * t19) - t27 * t1, t1 (t1 * t21 - t10 * t18) * r_i_i_C(1) + (-t1 * t18 - t10 * t21) * r_i_i_C(2); 0, t24 * t8 - t25 * t7, t7, -t26 * (t16 * t31 - t7 * t19) + t27 * t3, -t3 (-t8 * t18 - t3 * t21) * r_i_i_C(1) + (t3 * t18 - t8 * t21) * r_i_i_C(2); 1 (t24 * t20 + t25 * t23) * t15, -t30, t26 * (t17 * t22 - t19 * t30) - t27 * t11, t11 (t11 * t21 - t18 * t32) * r_i_i_C(1) + (-t11 * t18 - t21 * t32) * r_i_i_C(2);];
Ja_transl  = t2;
