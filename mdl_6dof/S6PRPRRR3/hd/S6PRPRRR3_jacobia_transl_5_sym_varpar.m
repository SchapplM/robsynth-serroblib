% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:48
% EndTime: 2019-02-26 19:54:48
% DurationCPUTime: 0.10s
% Computational Cost: add. (140->30), mult. (162->52), div. (0->0), fcn. (199->11), ass. (0->27)
t19 = pkin(12) + qJ(4);
t17 = qJ(5) + t19;
t13 = sin(t17);
t14 = cos(t17);
t21 = sin(pkin(6));
t22 = cos(pkin(11));
t30 = t21 * t22;
t20 = sin(pkin(11));
t25 = cos(qJ(2));
t23 = cos(pkin(6));
t24 = sin(qJ(2));
t28 = t23 * t24;
t8 = t20 * t25 + t22 * t28;
t35 = (-t8 * t13 - t14 * t30) * r_i_i_C(1) + (t13 * t30 - t8 * t14) * r_i_i_C(2);
t10 = -t20 * t28 + t22 * t25;
t31 = t20 * t21;
t34 = (-t10 * t13 + t14 * t31) * r_i_i_C(1) + (-t10 * t14 - t13 * t31) * r_i_i_C(2);
t29 = t21 * t24;
t33 = (-t13 * t29 + t23 * t14) * r_i_i_C(1) + (-t23 * t13 - t14 * t29) * r_i_i_C(2);
t32 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(3);
t27 = t23 * t25;
t16 = cos(t19);
t26 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + pkin(4) * t16 + cos(pkin(12)) * pkin(3) + pkin(2);
t15 = sin(t19);
t9 = t20 * t27 + t22 * t24;
t7 = t20 * t24 - t22 * t27;
t1 = [0, t32 * t10 - t26 * t9, t9 (-t10 * t15 + t16 * t31) * pkin(4) + t34, t34, 0; 0, -t26 * t7 + t32 * t8, t7 (-t15 * t8 - t16 * t30) * pkin(4) + t35, t35, 0; 1 (t32 * t24 + t26 * t25) * t21, -t21 * t25 (-t15 * t29 + t16 * t23) * pkin(4) + t33, t33, 0;];
Ja_transl  = t1;
