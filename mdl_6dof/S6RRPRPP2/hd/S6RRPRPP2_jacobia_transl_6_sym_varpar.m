% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPP2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:21
% EndTime: 2019-02-26 21:35:21
% DurationCPUTime: 0.11s
% Computational Cost: add. (154->29), mult. (186->37), div. (0->0), fcn. (216->8), ass. (0->24)
t21 = pkin(8) - r_i_i_C(3) - qJ(6);
t11 = qJ(2) + pkin(9);
t8 = sin(t11);
t30 = cos(qJ(2)) * pkin(2) + t21 * t8;
t9 = cos(t11);
t29 = t9 * pkin(3) + pkin(1) + t30;
t13 = sin(qJ(4));
t16 = cos(qJ(4));
t22 = pkin(4) + pkin(5) + r_i_i_C(1);
t23 = r_i_i_C(2) + qJ(5);
t28 = t23 * t13 + t22 * t16 + pkin(3);
t15 = sin(qJ(1));
t27 = t15 * t13;
t26 = t15 * t16;
t17 = cos(qJ(1));
t25 = t17 * t13;
t24 = t17 * t16;
t18 = -sin(qJ(2)) * pkin(2) + t21 * t9 - t28 * t8;
t12 = -qJ(3) - pkin(7);
t4 = t9 * t24 + t27;
t3 = t9 * t25 - t26;
t2 = t9 * t26 - t25;
t1 = t9 * t27 + t24;
t5 = [-t23 * t1 - t17 * t12 - t29 * t15 - t22 * t2, t18 * t17, t15, -t22 * t3 + t23 * t4, t3, -t17 * t8; -t15 * t12 + t29 * t17 + t22 * t4 + t23 * t3, t18 * t15, -t17, -t22 * t1 + t23 * t2, t1, -t15 * t8; 0, t28 * t9 + t30, 0 (-t22 * t13 + t23 * t16) * t8, t8 * t13, t9;];
Ja_transl  = t5;
