% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:46:35
% EndTime: 2019-02-26 21:46:35
% DurationCPUTime: 0.13s
% Computational Cost: add. (241->33), mult. (201->41), div. (0->0), fcn. (224->10), ass. (0->31)
t27 = sin(qJ(5));
t29 = cos(qJ(5));
t37 = r_i_i_C(3) + qJ(6);
t48 = pkin(5) + r_i_i_C(1);
t54 = t37 * t27 + t48 * t29;
t53 = pkin(9) + r_i_i_C(2);
t26 = qJ(2) + pkin(10);
t23 = qJ(4) + t26;
t20 = cos(t23);
t51 = t20 * t53;
t19 = sin(t23);
t50 = t20 * pkin(4) + t53 * t19;
t35 = pkin(3) * cos(t26) + cos(qJ(2)) * pkin(2);
t49 = pkin(1) + t35 + t50;
t28 = sin(qJ(1));
t46 = t28 * t51;
t41 = t28 * t27;
t40 = t28 * t29;
t30 = cos(qJ(1));
t39 = t30 * t27;
t38 = t30 * t29;
t36 = t30 * t51;
t33 = t54 * t20 + t50;
t32 = (-pkin(4) - t54) * t19;
t31 = -pkin(3) * sin(t26) - sin(qJ(2)) * pkin(2) + t32;
t25 = -pkin(8) - qJ(3) - pkin(7);
t4 = t20 * t38 + t41;
t3 = t20 * t39 - t40;
t2 = t20 * t40 - t39;
t1 = t20 * t41 + t38;
t5 = [-t37 * t1 - t48 * t2 - t30 * t25 - t49 * t28, t31 * t30 + t36, t28, t30 * t32 + t36, -t48 * t3 + t37 * t4, t3; -t28 * t25 + t37 * t3 + t49 * t30 + t48 * t4, t31 * t28 + t46, -t30, t28 * t32 + t46, -t48 * t1 + t37 * t2, t1; 0, t33 + t35, 0, t33 (-t48 * t27 + t37 * t29) * t19, t19 * t27;];
Ja_transl  = t5;
