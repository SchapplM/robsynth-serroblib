% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR3_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR3_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobia_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:33
% EndTime: 2019-02-26 21:29:33
% DurationCPUTime: 0.12s
% Computational Cost: add. (112->37), mult. (271->57), div. (0->0), fcn. (353->10), ass. (0->28)
t18 = sin(pkin(11));
t21 = cos(pkin(11));
t23 = sin(qJ(2));
t25 = cos(qJ(2));
t12 = t23 * t18 - t25 * t21;
t37 = t25 * pkin(2);
t22 = cos(pkin(6));
t36 = t22 * t25;
t19 = sin(pkin(6));
t24 = sin(qJ(1));
t34 = t24 * t19;
t26 = cos(qJ(1));
t32 = t26 * t19;
t31 = r_i_i_C(3) + qJ(4);
t29 = t25 * t18 + t23 * t21;
t10 = t29 * t22;
t3 = -t26 * t10 + t24 * t12;
t30 = t24 * t10 + t26 * t12;
t17 = sin(pkin(12));
t20 = cos(pkin(12));
t28 = t20 * r_i_i_C(1) - t17 * r_i_i_C(2) + pkin(3);
t27 = t12 * t22;
t16 = pkin(1) + t37;
t11 = t22 * t23 * pkin(2) + (-pkin(8) - qJ(3)) * t19;
t8 = t12 * t19;
t5 = t24 * t27 - t26 * t29;
t2 = -t24 * t29 - t26 * t27;
t1 = [(t17 * t32 + t3 * t20) * r_i_i_C(1) + (-t3 * t17 + t20 * t32) * r_i_i_C(2) + t3 * pkin(3) - t24 * t16 - t26 * t11 + t31 * t2, -t31 * t30 + (-t23 * t26 - t24 * t36) * pkin(2) + t28 * t5, t34, -t5, 0, 0; (t17 * t34 - t20 * t30) * r_i_i_C(1) + (t17 * t30 + t20 * t34) * r_i_i_C(2) - t30 * pkin(3) + t26 * t16 - t24 * t11 - t31 * t5, -t31 * t3 + (-t23 * t24 + t26 * t36) * pkin(2) + t28 * t2, -t32, -t2, 0, 0; 0, -t28 * t8 + (t29 * t31 + t37) * t19, t22, t8, 0, 0;];
Ja_transl  = t1;
