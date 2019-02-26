% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRP3_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP3_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:02:22
% EndTime: 2019-02-26 20:02:22
% DurationCPUTime: 0.10s
% Computational Cost: add. (84->23), mult. (222->43), div. (0->0), fcn. (279->10), ass. (0->24)
t15 = sin(pkin(6));
t19 = sin(qJ(3));
t30 = t15 * t19;
t21 = cos(qJ(3));
t29 = t15 * t21;
t18 = cos(pkin(6));
t20 = sin(qJ(2));
t28 = t18 * t20;
t22 = cos(qJ(2));
t27 = t18 * t22;
t26 = r_i_i_C(3) + qJ(4);
t13 = sin(pkin(11));
t16 = cos(pkin(11));
t25 = r_i_i_C(1) * t16 - r_i_i_C(2) * t13 + pkin(3);
t24 = t13 * r_i_i_C(1) + t16 * r_i_i_C(2) + pkin(8);
t23 = t26 * t19 + t25 * t21 + pkin(2);
t17 = cos(pkin(10));
t14 = sin(pkin(10));
t9 = -t18 * t21 + t20 * t30;
t8 = -t14 * t28 + t17 * t22;
t6 = t14 * t22 + t17 * t28;
t3 = -t14 * t29 + t19 * t8;
t1 = t17 * t29 + t19 * t6;
t2 = [0, t24 * t8 + t23 * (-t14 * t27 - t17 * t20) t26 * (t14 * t30 + t21 * t8) - t25 * t3, t3, 0, 0; 0, t24 * t6 + t23 * (-t14 * t20 + t17 * t27) t26 * (-t17 * t30 + t21 * t6) - t25 * t1, t1, 0, 0; 1 (t24 * t20 + t23 * t22) * t15, t26 * (t18 * t19 + t20 * t29) - t25 * t9, t9, 0, 0;];
Ja_transl  = t2;
