% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR10V2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10V2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobia_transl_5_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:31
% EndTime: 2019-04-11 14:56:31
% DurationCPUTime: 0.16s
% Computational Cost: add. (163->43), mult. (225->70), div. (0->0), fcn. (263->10), ass. (0->39)
t24 = qJ(2) + qJ(3);
t22 = cos(t24);
t21 = sin(t24);
t29 = cos(qJ(5));
t25 = sin(qJ(5));
t30 = cos(qJ(4));
t48 = t25 * t30;
t35 = t21 * t48 + t22 * t29;
t45 = t29 * t30;
t36 = -t21 * t45 + t22 * t25;
t54 = pkin(5) * t22 + r_i_i_C(1) * t36 + r_i_i_C(2) * t35;
t26 = sin(qJ(4));
t52 = r_i_i_C(3) * t26;
t16 = t21 * pkin(5);
t17 = t22 * pkin(3);
t51 = t21 * t25;
t50 = t21 * t29;
t31 = cos(qJ(1));
t49 = t21 * t31;
t28 = sin(qJ(1));
t47 = t28 * t26;
t46 = t28 * t30;
t44 = t31 * t26;
t43 = t31 * t30;
t42 = t54 * t28;
t41 = t54 * t31;
t23 = cos(qJ(2)) * pkin(2);
t40 = -t23 - pkin(1) - t17;
t39 = (t22 * t45 + t51) * r_i_i_C(1) + (-t22 * t48 + t50) * r_i_i_C(2) + t22 * t52 + t17 + t16;
t38 = (-pkin(3) - t52) * t21;
t37 = t29 * r_i_i_C(1) - t25 * r_i_i_C(2);
t32 = -sin(qJ(2)) * pkin(2) + t38;
t12 = t22 * t43 + t47;
t11 = t22 * t44 - t46;
t10 = t22 * t46 - t44;
t9 = -t22 * t47 - t43;
t2 = t12 * t29 + t25 * t49;
t1 = -t12 * t25 + t29 * t49;
t3 = [t9 * r_i_i_C(3) - t37 * t10 + ((-t25 * r_i_i_C(1) - t29 * r_i_i_C(2) - pkin(5)) * t21 + t40) * t28, t32 * t31 + t41, t31 * t38 + t41, t12 * r_i_i_C(3) - t37 * t11, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t11 * r_i_i_C(3) + (-t40 + t16) * t31, t32 * t28 + t42, t28 * t38 + t42, t10 * r_i_i_C(3) + t37 * t9 (-t10 * t25 + t28 * t50) * r_i_i_C(1) + (-t10 * t29 - t28 * t51) * r_i_i_C(2), 0; 0, t23 + t39, t39 (r_i_i_C(3) * t30 - t37 * t26) * t21, -t35 * r_i_i_C(1) + t36 * r_i_i_C(2), 0;];
Ja_transl  = t3;
