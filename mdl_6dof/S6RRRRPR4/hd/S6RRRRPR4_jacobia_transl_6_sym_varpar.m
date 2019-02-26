% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:32:17
% EndTime: 2019-02-26 22:32:17
% DurationCPUTime: 0.16s
% Computational Cost: add. (261->42), mult. (182->54), div. (0->0), fcn. (197->12), ass. (0->39)
t29 = qJ(4) + pkin(11);
t15 = pkin(5) * cos(t29) + cos(qJ(4)) * pkin(4);
t12 = pkin(3) + t15;
t30 = qJ(2) + qJ(3);
t24 = sin(t30);
t25 = cos(t30);
t28 = -pkin(10) - qJ(5) - pkin(9);
t56 = t25 * t12 + (r_i_i_C(3) - t28) * t24;
t27 = cos(qJ(2)) * pkin(2);
t55 = pkin(1) + t27 + t56;
t23 = qJ(6) + t29;
t20 = cos(t23);
t33 = cos(qJ(1));
t44 = t33 * t20;
t19 = sin(t23);
t32 = sin(qJ(1));
t47 = t32 * t19;
t5 = t25 * t47 + t44;
t45 = t33 * t19;
t46 = t32 * t20;
t6 = -t25 * t46 + t45;
t54 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t25 * t45 + t46;
t8 = t25 * t44 + t47;
t53 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t52 = r_i_i_C(1) * t20;
t51 = r_i_i_C(2) * t19;
t49 = t25 * t32;
t48 = t25 * t33;
t40 = t24 * t51;
t43 = r_i_i_C(3) * t49 + t32 * t40;
t42 = r_i_i_C(3) * t48 + t33 * t40;
t14 = pkin(5) * sin(t29) + sin(qJ(4)) * pkin(4);
t41 = t14 + pkin(8) + pkin(7);
t38 = -r_i_i_C(1) * t19 - r_i_i_C(2) * t20;
t37 = -t25 * t28 + (-t12 - t52) * t24;
t36 = (-t51 + t52) * t25 + t56;
t35 = -sin(qJ(2)) * pkin(2) + t37;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t55 * t32 + t41 * t33, t35 * t33 + t42, t37 * t33 + t42, -t14 * t48 + t32 * t15 + t53, t33 * t24, t53; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t41 * t32 + t55 * t33, t35 * t32 + t43, t37 * t32 + t43, -t14 * t49 - t33 * t15 + t54, t32 * t24, t54; 0, t27 + t36, t36 (-t14 + t38) * t24, -t25, t38 * t24;];
Ja_transl  = t1;
