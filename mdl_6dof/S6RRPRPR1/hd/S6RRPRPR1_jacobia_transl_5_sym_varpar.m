% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:40
% EndTime: 2019-02-26 21:37:40
% DurationCPUTime: 0.12s
% Computational Cost: add. (147->24), mult. (108->25), div. (0->0), fcn. (117->10), ass. (0->21)
t43 = r_i_i_C(3) + qJ(5);
t19 = qJ(2) + pkin(10);
t16 = qJ(4) + t19;
t14 = sin(t16);
t15 = cos(t16);
t21 = cos(pkin(11));
t29 = -r_i_i_C(1) * t21 - pkin(4);
t20 = sin(pkin(11));
t37 = r_i_i_C(2) * t20;
t24 = t43 * t14 + (-t29 - t37) * t15;
t42 = t24 + pkin(3) * cos(t19) + cos(qJ(2)) * pkin(2);
t41 = t14 * t37 + t43 * t15;
t39 = pkin(1) + t42;
t22 = sin(qJ(1));
t32 = t41 * t22;
t23 = cos(qJ(1));
t31 = t41 * t23;
t28 = t29 * t14;
t26 = t20 * r_i_i_C(1) + t21 * r_i_i_C(2) + pkin(7) + pkin(8) + qJ(3);
t25 = -pkin(3) * sin(t19) - sin(qJ(2)) * pkin(2) + t28;
t1 = [-t39 * t22 + t26 * t23, t25 * t23 + t31, t22, t23 * t28 + t31, t23 * t14, 0; t26 * t22 + t39 * t23, t25 * t22 + t32, -t23, t22 * t28 + t32, t22 * t14, 0; 0, t42, 0, t24, -t15, 0;];
Ja_transl  = t1;
