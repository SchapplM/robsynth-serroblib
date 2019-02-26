% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

function Ja_transl = S6PRPPRR1_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR1_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobia_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:35
% EndTime: 2019-02-26 19:44:35
% DurationCPUTime: 0.11s
% Computational Cost: add. (60->20), mult. (155->34), div. (0->0), fcn. (203->10), ass. (0->19)
t13 = sin(pkin(11));
t17 = cos(pkin(11));
t20 = sin(qJ(2));
t21 = cos(qJ(2));
t9 = t20 * t13 - t17 * t21;
t19 = cos(pkin(6));
t26 = t19 * t21;
t24 = r_i_i_C(3) + qJ(4);
t23 = t13 * t21 + t20 * t17;
t22 = r_i_i_C(1) * cos(pkin(12)) - r_i_i_C(2) * sin(pkin(12)) + pkin(3);
t18 = cos(pkin(10));
t15 = sin(pkin(6));
t14 = sin(pkin(10));
t8 = t23 * t19;
t7 = t9 * t19;
t5 = t9 * t15;
t4 = t14 * t7 - t18 * t23;
t2 = -t14 * t23 - t18 * t7;
t1 = [0, -t24 * (t14 * t8 + t18 * t9) + (-t14 * t26 - t18 * t20) * pkin(2) + t22 * t4, t14 * t15, -t4, 0, 0; 0, -t24 * (t14 * t9 - t18 * t8) + (-t14 * t20 + t18 * t26) * pkin(2) + t22 * t2, -t18 * t15, -t2, 0, 0; 1, -t22 * t5 + (pkin(2) * t21 + t23 * t24) * t15, t19, t5, 0, 0;];
Ja_transl  = t1;
