% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:10
% EndTime: 2019-02-26 19:56:10
% DurationCPUTime: 0.10s
% Computational Cost: add. (102->28), mult. (179->51), div. (0->0), fcn. (220->10), ass. (0->25)
t14 = qJ(4) + qJ(5);
t12 = sin(t14);
t13 = cos(t14);
t15 = sin(pkin(11));
t16 = sin(pkin(6));
t31 = t15 * t16;
t17 = cos(pkin(11));
t20 = sin(qJ(2));
t18 = cos(pkin(6));
t22 = cos(qJ(2));
t26 = t18 * t22;
t9 = t15 * t26 + t17 * t20;
t34 = (-t12 * t31 + t13 * t9) * r_i_i_C(1) + (-t12 * t9 - t13 * t31) * r_i_i_C(2);
t30 = t16 * t17;
t7 = t15 * t20 - t17 * t26;
t33 = (t12 * t30 + t13 * t7) * r_i_i_C(1) + (-t12 * t7 + t13 * t30) * r_i_i_C(2);
t28 = t16 * t22;
t32 = (-t12 * t18 - t13 * t28) * r_i_i_C(1) + (t12 * t28 - t13 * t18) * r_i_i_C(2);
t19 = sin(qJ(4));
t29 = t16 * t19;
t27 = t18 * t20;
t25 = pkin(2) + r_i_i_C(3) + pkin(9) + pkin(8);
t24 = pkin(4) * t19 + r_i_i_C(1) * t12 + r_i_i_C(2) * t13 + qJ(3);
t21 = cos(qJ(4));
t1 = [0, -t25 * t9 + t24 * (-t15 * t27 + t17 * t22) t9 (-t15 * t29 + t21 * t9) * pkin(4) + t34, t34, 0; 0, -t25 * t7 + t24 * (t15 * t22 + t17 * t27) t7 (t17 * t29 + t21 * t7) * pkin(4) + t33, t33, 0; 1 (t24 * t20 + t25 * t22) * t16, -t28 (-t18 * t19 - t21 * t28) * pkin(4) + t32, t32, 0;];
Ja_transl  = t1;
