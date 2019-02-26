% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR13_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR13_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:49
% EndTime: 2019-02-26 21:44:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (146->37), mult. (363->57), div. (0->0), fcn. (461->10), ass. (0->30)
t17 = sin(pkin(6));
t21 = sin(qJ(1));
t36 = t17 * t21;
t23 = cos(qJ(2));
t35 = t17 * t23;
t24 = cos(qJ(1));
t34 = t17 * t24;
t33 = r_i_i_C(3) + qJ(5);
t32 = cos(pkin(6));
t31 = t17 * (pkin(3) + pkin(8));
t30 = t21 * t32;
t29 = t24 * t32;
t16 = sin(pkin(11));
t18 = cos(pkin(11));
t28 = r_i_i_C(1) * t18 - r_i_i_C(2) * t16 + pkin(4);
t20 = sin(qJ(2));
t10 = t21 * t20 - t23 * t29;
t19 = sin(qJ(4));
t22 = cos(qJ(4));
t27 = -t10 * t19 + t22 * t34;
t3 = t10 * t22 + t19 * t34;
t26 = t16 * r_i_i_C(1) + t18 * r_i_i_C(2) + pkin(2) + pkin(9);
t25 = t28 * t19 - t33 * t22 + qJ(3);
t13 = -t20 * t30 + t24 * t23;
t12 = t24 * t20 + t23 * t30;
t11 = t20 * t29 + t21 * t23;
t8 = t32 * t19 + t22 * t35;
t2 = t12 * t19 + t22 * t36;
t1 = -t12 * t22 + t19 * t36;
t4 = [-t21 * pkin(1) - t10 * qJ(3) - t26 * t11 + t24 * t31 + t28 * t27 + t33 * t3, -t26 * t12 + t25 * t13, t12, -t28 * t1 + t33 * t2, t1, 0; t24 * pkin(1) + t12 * qJ(3) + t33 * t1 + t26 * t13 + t28 * t2 + t21 * t31, -t26 * t10 + t25 * t11, t10, -t33 * t27 + t28 * t3, -t3, 0; 0 (t25 * t20 + t26 * t23) * t17, -t35, t33 * (-t19 * t35 + t32 * t22) - t28 * t8, t8, 0;];
Ja_transl  = t4;
