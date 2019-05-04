% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR14V3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14V3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobia_transl_6_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:09
% EndTime: 2019-04-12 15:12:09
% DurationCPUTime: 0.18s
% Computational Cost: add. (117->53), mult. (330->105), div. (0->0), fcn. (419->10), ass. (0->42)
t23 = sin(qJ(4));
t28 = cos(qJ(4));
t30 = cos(qJ(1));
t37 = t30 * t28;
t25 = sin(qJ(1));
t29 = cos(qJ(2));
t43 = t25 * t29;
t12 = t23 * t43 + t37;
t21 = sin(qJ(6));
t26 = cos(qJ(6));
t39 = t30 * t23;
t13 = t28 * t43 - t39;
t22 = sin(qJ(5));
t27 = cos(qJ(5));
t24 = sin(qJ(2));
t44 = t25 * t24;
t4 = t13 * t27 + t22 * t44;
t54 = -t12 * t26 + t4 * t21;
t53 = -t12 * t21 - t4 * t26;
t41 = t29 * t22;
t45 = t24 * t28;
t11 = t27 * t45 - t41;
t40 = t29 * t27;
t33 = -t22 * t45 - t40;
t34 = t26 * r_i_i_C(1) - t21 * r_i_i_C(2);
t47 = t23 * t24;
t52 = t33 * r_i_i_C(3) + t29 * qJ(3) + (-t21 * r_i_i_C(1) - t26 * r_i_i_C(2)) * t47 - t34 * t11;
t51 = t22 * r_i_i_C(3);
t48 = t21 * t27;
t46 = t23 * t29;
t42 = t26 * t27;
t38 = t30 * t24;
t36 = t24 * qJ(3);
t35 = -t13 * t22 + t27 * t44;
t16 = t25 * t23 + t29 * t37;
t15 = -t25 * t28 + t29 * t39;
t14 = t24 * t22 + t28 * t40;
t7 = t16 * t27 + t22 * t38;
t6 = t16 * t22 - t27 * t38;
t2 = t15 * t21 + t7 * t26;
t1 = t15 * t26 - t7 * t21;
t3 = [t53 * r_i_i_C(1) + t54 * r_i_i_C(2) + t35 * r_i_i_C(3) - t25 * t36, t52 * t30, t38 (-t15 * t42 + t16 * t21) * r_i_i_C(1) + (t15 * t48 + t16 * t26) * r_i_i_C(2) - t15 * t51, t7 * r_i_i_C(3) - t34 * t6, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t6 * r_i_i_C(3) + t30 * t36, t52 * t25, t44 (-t12 * t42 + t13 * t21) * r_i_i_C(1) + (t12 * t48 + t13 * t26) * r_i_i_C(2) - t12 * t51, t4 * r_i_i_C(3) + t34 * t35, -t54 * r_i_i_C(1) + t53 * r_i_i_C(2); 0 (t14 * t26 + t21 * t46) * r_i_i_C(1) + (-t14 * t21 + t26 * t46) * r_i_i_C(2) + (-t24 * t27 + t28 * t41) * r_i_i_C(3) + t36, -t29 ((t21 * t28 - t23 * t42) * r_i_i_C(1) + (t23 * t48 + t26 * t28) * r_i_i_C(2) - t23 * t51) * t24, t11 * r_i_i_C(3) + t34 * t33 (-t11 * t21 + t26 * t47) * r_i_i_C(1) + (-t11 * t26 - t21 * t47) * r_i_i_C(2);];
Ja_transl  = t3;
