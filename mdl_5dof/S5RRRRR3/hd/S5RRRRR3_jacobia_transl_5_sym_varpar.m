% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
%
% Output:
% Ja_transl [3x5]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-31 10:31
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobia_transl_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobia_transl_5_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-31 10:31:13
% EndTime: 2019-05-31 10:31:13
% DurationCPUTime: 0.16s
% Computational Cost: add. (165->35), mult. (163->52), div. (0->0), fcn. (175->10), ass. (0->36)
t54 = pkin(5) + r_i_i_C(3);
t29 = cos(qJ(4));
t19 = pkin(3) * t29 + pkin(2);
t25 = qJ(2) + qJ(3);
t21 = sin(t25);
t23 = cos(t25);
t53 = t23 * t19 + t54 * t21;
t48 = cos(qJ(2)) * pkin(1);
t52 = t48 + t53;
t24 = qJ(4) + qJ(5);
t22 = cos(t24);
t31 = cos(qJ(1));
t20 = sin(t24);
t28 = sin(qJ(1));
t41 = t28 * t20;
t5 = t22 * t31 + t23 * t41;
t40 = t28 * t22;
t6 = t20 * t31 - t23 * t40;
t50 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t44 = t23 * t31;
t7 = -t20 * t44 + t40;
t8 = t22 * t44 + t41;
t49 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t47 = r_i_i_C(1) * t22;
t46 = r_i_i_C(2) * t20;
t26 = sin(qJ(4));
t43 = t26 * t28;
t42 = t26 * t31;
t39 = t21 * t46;
t38 = (t23 * t54 + t39) * t28;
t37 = t31 * t39 + t54 * t44;
t36 = -r_i_i_C(1) * t20 - r_i_i_C(2) * t22;
t35 = (-t19 - t47) * t21;
t33 = (-t46 + t47) * t23 + t53;
t32 = -sin(qJ(2)) * pkin(1) + t35;
t1 = [pkin(3) * t42 + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t52 * t28, t32 * t31 + t37, t31 * t35 + t37, (-t23 * t42 + t28 * t29) * pkin(3) + t49, t49; pkin(3) * t43 + t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t52 * t31, t32 * t28 + t38, t28 * t35 + t38, (-t23 * t43 - t29 * t31) * pkin(3) + t50, t50; 0, t33 + t48, t33, (-pkin(3) * t26 + t36) * t21, t36 * t21;];
Ja_transl  = t1;
