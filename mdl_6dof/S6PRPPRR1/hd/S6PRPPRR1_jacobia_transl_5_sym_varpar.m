% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPPRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:35
% EndTime: 2019-02-26 19:44:35
% DurationCPUTime: 0.12s
% Computational Cost: add. (107->31), mult. (220->54), div. (0->0), fcn. (288->11), ass. (0->27)
t19 = sin(pkin(11));
t22 = cos(pkin(11));
t26 = sin(qJ(2));
t27 = cos(qJ(2));
t11 = t26 * t19 - t27 * t22;
t36 = r_i_i_C(3) + pkin(8) + qJ(4);
t20 = sin(pkin(10));
t21 = sin(pkin(6));
t35 = t20 * t21;
t23 = cos(pkin(10));
t34 = t23 * t21;
t24 = cos(pkin(6));
t33 = t24 * t27;
t29 = t27 * t19 + t26 * t22;
t10 = t29 * t24;
t3 = t23 * t10 - t20 * t11;
t30 = t20 * t10 + t23 * t11;
t18 = pkin(12) + qJ(5);
t16 = sin(t18);
t17 = cos(t18);
t28 = t17 * r_i_i_C(1) - t16 * r_i_i_C(2) + cos(pkin(12)) * pkin(4) + pkin(3);
t9 = t11 * t24;
t8 = t29 * t21;
t7 = t11 * t21;
t5 = t20 * t9 - t23 * t29;
t2 = -t20 * t29 - t23 * t9;
t1 = [0, -t36 * t30 + (-t20 * t33 - t23 * t26) * pkin(2) + t28 * t5, t35, -t5 (t16 * t30 + t17 * t35) * r_i_i_C(1) + (-t16 * t35 + t17 * t30) * r_i_i_C(2), 0; 0, t36 * t3 + (-t20 * t26 + t23 * t33) * pkin(2) + t28 * t2, -t34, -t2 (-t3 * t16 - t17 * t34) * r_i_i_C(1) + (t16 * t34 - t3 * t17) * r_i_i_C(2), 0; 1, t21 * t27 * pkin(2) - t28 * t7 + t36 * t8, t24, t7 (-t8 * t16 + t24 * t17) * r_i_i_C(1) + (-t24 * t16 - t8 * t17) * r_i_i_C(2), 0;];
Ja_transl  = t1;
