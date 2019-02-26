% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:09:30
% EndTime: 2019-02-26 22:09:30
% DurationCPUTime: 0.15s
% Computational Cost: add. (201->39), mult. (157->46), div. (0->0), fcn. (170->10), ass. (0->32)
t22 = qJ(2) + qJ(3);
t18 = pkin(10) + t22;
t14 = sin(t18);
t15 = cos(t18);
t24 = sin(qJ(5));
t43 = r_i_i_C(2) * t24;
t50 = r_i_i_C(3) * t15 + t14 * t43;
t26 = cos(qJ(5));
t17 = pkin(5) * t26 + pkin(4);
t23 = -qJ(6) - pkin(9);
t49 = (r_i_i_C(3) - t23) * t14 + t15 * t17 + pkin(3) * cos(t22);
t48 = pkin(5) + r_i_i_C(1);
t44 = r_i_i_C(1) * t26;
t28 = (-t17 - t44) * t14 - t15 * t23 - pkin(3) * sin(t22);
t20 = cos(qJ(2)) * pkin(2);
t46 = pkin(1) + t20 + t49;
t25 = sin(qJ(1));
t40 = t50 * t25;
t27 = cos(qJ(1));
t39 = t50 * t27;
t38 = t24 * t27;
t37 = t25 * t24;
t36 = t25 * t26;
t35 = t26 * t27;
t33 = pkin(5) * t24 + pkin(7) + pkin(8) + qJ(4);
t3 = -t15 * t38 + t36;
t1 = t15 * t37 + t35;
t30 = -sin(qJ(2)) * pkin(2) + t28;
t29 = (-t43 + t44) * t15 + t49;
t4 = t15 * t35 + t37;
t2 = -t15 * t36 + t38;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t46 * t25 + t33 * t27, t30 * t27 + t39, t28 * t27 + t39, t25, -t4 * r_i_i_C(2) + t48 * t3, t27 * t14; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t33 * t25 + t46 * t27, t30 * t25 + t40, t28 * t25 + t40, -t27, t2 * r_i_i_C(2) - t48 * t1, t25 * t14; 0, t20 + t29, t29, 0 (-r_i_i_C(2) * t26 - t48 * t24) * t14, -t15;];
Ja_transl  = t5;
