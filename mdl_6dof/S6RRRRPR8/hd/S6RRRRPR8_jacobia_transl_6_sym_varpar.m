% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:42
% EndTime: 2019-02-26 22:34:43
% DurationCPUTime: 0.21s
% Computational Cost: add. (312->49), mult. (388->71), div. (0->0), fcn. (471->10), ass. (0->40)
t54 = pkin(4) + pkin(5);
t32 = cos(qJ(3));
t23 = t32 * pkin(3) + pkin(2);
t33 = cos(qJ(2));
t29 = sin(qJ(2));
t48 = -r_i_i_C(3) - pkin(10) + pkin(9) + pkin(8);
t42 = t48 * t29;
t62 = t33 * t23 + pkin(1) + t42;
t26 = qJ(3) + qJ(4);
t24 = sin(t26);
t25 = cos(t26);
t34 = cos(qJ(1));
t49 = t34 * t25;
t30 = sin(qJ(1));
t51 = t30 * t33;
t16 = t24 * t51 + t49;
t50 = t34 * t24;
t17 = t25 * t51 - t50;
t27 = sin(qJ(6));
t31 = cos(qJ(6));
t9 = t16 * t31;
t61 = -(t17 * t27 - t9) * r_i_i_C(1) - (t16 * t27 + t17 * t31) * r_i_i_C(2);
t60 = (-(t24 * t27 + t25 * t31) * r_i_i_C(2) - (-t24 * t31 + t25 * t27) * r_i_i_C(1)) * t29;
t18 = -t30 * t25 + t33 * t50;
t19 = t30 * t24 + t33 * t49;
t5 = t18 * t31 - t19 * t27;
t6 = t18 * t27 + t19 * t31;
t58 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t37 = t31 * r_i_i_C(1) - t27 * r_i_i_C(2) + t54;
t43 = -t27 * r_i_i_C(1) - qJ(5);
t56 = (t31 * r_i_i_C(2) - t43) * t24 + t37 * t25 + t23;
t28 = sin(qJ(3));
t53 = t28 * pkin(3);
t47 = t29 * t25 * qJ(5) - t60;
t46 = t54 * t24;
t45 = pkin(7) + t53;
t39 = t17 * qJ(5) - t54 * t16 - t61;
t38 = t19 * qJ(5) - t54 * t18 - t58;
t36 = -t56 * t29 + t48 * t33;
t1 = [-t9 * r_i_i_C(2) + t43 * t16 - t37 * t17 - t62 * t30 + t45 * t34, t36 * t34 (-t28 * t33 * t34 + t30 * t32) * pkin(3) + t38, t38, t18, t58; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t18 * qJ(5) + t54 * t19 + t45 * t30 + t62 * t34, t36 * t30 (-t28 * t51 - t32 * t34) * pkin(3) + t39, t39, t16, t61; 0, t56 * t33 + t42 (-t46 - t53) * t29 + t47, -t29 * t46 + t47, t29 * t24, t60;];
Ja_transl  = t1;
