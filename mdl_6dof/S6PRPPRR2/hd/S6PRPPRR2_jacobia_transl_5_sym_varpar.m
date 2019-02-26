% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPPRR2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:12
% EndTime: 2019-02-26 19:45:12
% DurationCPUTime: 0.12s
% Computational Cost: add. (93->29), mult. (240->52), div. (0->0), fcn. (315->10), ass. (0->23)
t12 = sin(pkin(11));
t15 = cos(pkin(11));
t19 = sin(qJ(2));
t21 = cos(qJ(2));
t9 = t19 * t12 - t21 * t15;
t13 = sin(pkin(10));
t14 = sin(pkin(6));
t31 = t13 * t14;
t16 = cos(pkin(10));
t30 = t16 * t14;
t17 = cos(pkin(6));
t29 = t17 * t21;
t26 = pkin(3) + pkin(8) + r_i_i_C(3);
t25 = t21 * t12 + t19 * t15;
t18 = sin(qJ(5));
t20 = cos(qJ(5));
t24 = t18 * r_i_i_C(1) + t20 * r_i_i_C(2) + qJ(4);
t23 = t25 * t17;
t22 = t9 * t17;
t7 = t9 * t14;
t4 = t13 * t22 - t16 * t25;
t2 = -t13 * t25 - t16 * t22;
t1 = [0 (-t13 * t29 - t16 * t19) * pkin(2) + t26 * t4 - t24 * (t13 * t23 + t16 * t9) t31, -t4 (-t18 * t31 - t4 * t20) * r_i_i_C(1) + (t4 * t18 - t20 * t31) * r_i_i_C(2), 0; 0 (-t13 * t19 + t16 * t29) * pkin(2) + t26 * t2 - t24 * (t13 * t9 - t16 * t23) -t30, -t2 (t18 * t30 - t2 * t20) * r_i_i_C(1) + (t2 * t18 + t20 * t30) * r_i_i_C(2), 0; 1, -t26 * t7 + (pkin(2) * t21 + t24 * t25) * t14, t17, t7 (-t17 * t18 + t7 * t20) * r_i_i_C(1) + (-t17 * t20 - t7 * t18) * r_i_i_C(2), 0;];
Ja_transl  = t1;
