% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPPRR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:49
% EndTime: 2019-02-26 19:45:50
% DurationCPUTime: 0.11s
% Computational Cost: add. (90->35), mult. (222->59), div. (0->0), fcn. (291->10), ass. (0->25)
t32 = pkin(2) + pkin(3);
t31 = pkin(8) + r_i_i_C(3);
t17 = sin(pkin(10));
t18 = sin(pkin(6));
t30 = t17 * t18;
t20 = cos(pkin(10));
t29 = t20 * t18;
t21 = cos(pkin(6));
t23 = sin(qJ(2));
t28 = t21 * t23;
t25 = cos(qJ(2));
t27 = t21 * t25;
t11 = t17 * t23 - t20 * t27;
t12 = t17 * t25 + t20 * t28;
t16 = sin(pkin(11));
t19 = cos(pkin(11));
t3 = t11 * t16 + t12 * t19;
t13 = t17 * t27 + t20 * t23;
t14 = -t17 * t28 + t20 * t25;
t6 = t13 * t16 + t14 * t19;
t22 = sin(qJ(5));
t24 = cos(qJ(5));
t26 = t24 * r_i_i_C(1) - t22 * r_i_i_C(2) + pkin(4);
t10 = (-t16 * t25 + t19 * t23) * t18;
t1 = [0, t14 * qJ(3) - t31 * t6 - t32 * t13 + t26 * (-t13 * t19 + t14 * t16) t13, -t30 (-t6 * t22 - t24 * t30) * r_i_i_C(1) + (t22 * t30 - t6 * t24) * r_i_i_C(2), 0; 0, t12 * qJ(3) - t32 * t11 - t31 * t3 + t26 * (-t11 * t19 + t12 * t16) t11, t29 (-t3 * t22 + t24 * t29) * r_i_i_C(1) + (-t22 * t29 - t3 * t24) * r_i_i_C(2), 0; 1, -t31 * t10 + (t26 * (t16 * t23 + t19 * t25) + qJ(3) * t23 + t32 * t25) * t18, -t18 * t25, -t21 (-t10 * t22 - t21 * t24) * r_i_i_C(1) + (-t10 * t24 + t21 * t22) * r_i_i_C(2), 0;];
Ja_transl  = t1;
