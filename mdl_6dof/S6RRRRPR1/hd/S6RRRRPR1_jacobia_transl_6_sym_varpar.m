% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR1
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
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:46
% EndTime: 2019-02-26 22:30:46
% DurationCPUTime: 0.16s
% Computational Cost: add. (303->37), mult. (175->45), div. (0->0), fcn. (185->12), ass. (0->34)
t28 = qJ(2) + qJ(3);
t26 = qJ(4) + t28;
t22 = pkin(11) + t26;
t18 = sin(t22);
t19 = cos(t22);
t29 = sin(qJ(6));
t49 = r_i_i_C(2) * t29;
t57 = pkin(10) + r_i_i_C(3);
t58 = t18 * t49 + t19 * t57;
t31 = cos(qJ(6));
t50 = r_i_i_C(1) * t31;
t35 = (-pkin(5) - t50) * t18 - pkin(4) * sin(t26);
t55 = t57 * t18 + t19 * pkin(5) + pkin(4) * cos(t26);
t36 = -pkin(3) * sin(t28) + t35;
t21 = pkin(3) * cos(t28);
t27 = cos(qJ(2)) * pkin(2);
t53 = pkin(1) + t21 + t27 + t55;
t32 = cos(qJ(1));
t46 = t29 * t32;
t30 = sin(qJ(1));
t45 = t30 * t29;
t44 = t30 * t31;
t43 = t31 * t32;
t41 = t58 * t30;
t40 = t58 * t32;
t37 = -sin(qJ(2)) * pkin(2) + t36;
t34 = (-t49 + t50) * t19 + t55;
t33 = t21 + t34;
t25 = -qJ(5) - pkin(9) - pkin(8) - pkin(7);
t4 = t19 * t43 + t45;
t3 = -t19 * t46 + t44;
t2 = -t19 * t44 + t46;
t1 = t19 * t45 + t43;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t25 * t32 - t53 * t30, t37 * t32 + t40, t36 * t32 + t40, t35 * t32 + t40, t30, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t30 * t25 + t53 * t32, t37 * t30 + t41, t36 * t30 + t41, t35 * t30 + t41, -t32, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, t27 + t33, t33, t34, 0 (-r_i_i_C(1) * t29 - r_i_i_C(2) * t31) * t18;];
Ja_transl  = t5;
