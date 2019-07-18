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
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
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
% StartTime: 2019-07-18 17:19:38
% EndTime: 2019-07-18 17:19:38
% DurationCPUTime: 0.16s
% Computational Cost: add. (165->35), mult. (163->52), div. (0->0), fcn. (175->10), ass. (0->37)
t55 = pkin(5) + r_i_i_C(3);
t29 = cos(qJ(4));
t19 = pkin(3) * t29 + pkin(2);
t25 = qJ(2) + qJ(3);
t21 = sin(t25);
t23 = cos(t25);
t54 = t23 * t19 + t55 * t21;
t49 = cos(qJ(2)) * pkin(1);
t53 = t49 + t54;
t24 = qJ(4) + qJ(5);
t22 = cos(t24);
t31 = cos(qJ(1));
t20 = sin(t24);
t28 = sin(qJ(1));
t42 = t28 * t20;
t5 = t22 * t31 + t23 * t42;
t40 = t31 * t20;
t41 = t28 * t22;
t6 = -t23 * t41 + t40;
t51 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t23 * t40 + t41;
t45 = t23 * t31;
t8 = t22 * t45 + t42;
t50 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t48 = r_i_i_C(1) * t22;
t47 = r_i_i_C(2) * t20;
t26 = sin(qJ(4));
t44 = t26 * t28;
t43 = t26 * t31;
t39 = t21 * t47;
t38 = (t23 * t55 + t39) * t28;
t37 = t31 * t39 + t55 * t45;
t36 = -r_i_i_C(1) * t20 - r_i_i_C(2) * t22;
t35 = (-t19 - t48) * t21;
t33 = (-t47 + t48) * t23 + t54;
t32 = -sin(qJ(2)) * pkin(1) + t35;
t1 = [pkin(3) * t43 + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t28 * t53, t32 * t31 + t37, t31 * t35 + t37, (-t23 * t43 + t28 * t29) * pkin(3) + t50, t50; pkin(3) * t44 + t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t31 * t53, t32 * t28 + t38, t28 * t35 + t38, (-t23 * t44 - t29 * t31) * pkin(3) + t51, t51; 0, t33 + t49, t33, (-pkin(3) * t26 + t36) * t21, t36 * t21;];
Ja_transl  = t1;
