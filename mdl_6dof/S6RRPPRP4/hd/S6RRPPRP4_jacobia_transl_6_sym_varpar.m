% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRP4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:46
% EndTime: 2019-02-26 21:26:47
% DurationCPUTime: 0.16s
% Computational Cost: add. (143->38), mult. (357->56), div. (0->0), fcn. (445->8), ass. (0->32)
t23 = sin(pkin(9));
t24 = cos(pkin(9));
t28 = cos(qJ(5));
t25 = sin(qJ(5));
t34 = t23 * t25 + t24 * t28;
t42 = r_i_i_C(3) + qJ(6);
t47 = t24 * t25;
t48 = pkin(5) + r_i_i_C(1);
t49 = pkin(3) + pkin(4);
t55 = -qJ(4) * t23 - t49 * t24 - t48 * t34 - pkin(2) - t42 * (-t23 * t28 + t47);
t29 = cos(qJ(2));
t26 = sin(qJ(2));
t41 = pkin(8) + r_i_i_C(2) - qJ(3);
t38 = t41 * t26;
t54 = -t29 * pkin(2) - pkin(1) + t38;
t51 = t55 * t26 - t41 * t29;
t46 = t26 * t23;
t27 = sin(qJ(1));
t45 = t27 * t29;
t30 = cos(qJ(1));
t44 = t30 * t23;
t43 = t30 * t24;
t17 = t23 * t45 + t43;
t18 = t24 * t45 - t44;
t37 = t17 * t28 - t18 * t25;
t36 = t17 * t25 + t18 * t28;
t20 = t27 * t23 + t29 * t43;
t19 = -t27 * t24 + t29 * t44;
t13 = t26 * t47 - t28 * t46;
t6 = t19 * t25 + t20 * t28;
t5 = -t19 * t28 + t20 * t25;
t1 = [t30 * pkin(7) - t17 * qJ(4) - t49 * t18 + t54 * t27 - t48 * t36 + t42 * t37, t51 * t30, t30 * t26, t19, t42 * t6 - t48 * t5, t5; t27 * pkin(7) + t19 * qJ(4) + t49 * t20 - t54 * t30 + t42 * t5 + t48 * t6, t51 * t27, t27 * t26, t17, t42 * t36 + t37 * t48, -t37; 0, -t55 * t29 - t38, -t29, t46, t42 * t34 * t26 - t48 * t13, t13;];
Ja_transl  = t1;
