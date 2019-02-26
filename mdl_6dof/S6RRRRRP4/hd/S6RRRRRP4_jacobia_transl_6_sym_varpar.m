% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:41:10
% EndTime: 2019-02-26 22:41:10
% DurationCPUTime: 0.18s
% Computational Cost: add. (297->39), mult. (260->50), div. (0->0), fcn. (289->10), ass. (0->38)
t30 = qJ(4) + qJ(5);
t25 = sin(t30);
t27 = cos(t30);
t48 = r_i_i_C(3) + qJ(6);
t58 = pkin(5) + r_i_i_C(1);
t62 = t48 * t25 + t58 * t27;
t35 = cos(qJ(4));
t23 = pkin(4) * t35 + pkin(3);
t31 = qJ(2) + qJ(3);
t26 = sin(t31);
t28 = cos(t31);
t37 = -pkin(10) - pkin(9);
t60 = t28 * t23 + (r_i_i_C(2) - t37) * t26;
t29 = cos(qJ(2)) * pkin(2);
t59 = pkin(1) + t29 + t60;
t32 = sin(qJ(4));
t57 = pkin(4) * t32;
t34 = sin(qJ(1));
t52 = t28 * t34;
t36 = cos(qJ(1));
t51 = t28 * t36;
t50 = t34 * t27;
t49 = t36 * t25;
t47 = t48 * t26 * t27;
t46 = t58 * t25;
t45 = pkin(8) + pkin(7) + t57;
t7 = t25 * t52 + t27 * t36;
t8 = t28 * t50 - t49;
t43 = t48 * t8 - t58 * t7;
t10 = t25 * t34 + t27 * t51;
t9 = t28 * t49 - t50;
t42 = t48 * t10 - t58 * t9;
t41 = t62 * t28 + t60;
t40 = -t28 * t37 + (-t23 - t62) * t26;
t39 = -sin(qJ(2)) * pkin(2) + t40;
t19 = r_i_i_C(2) * t51;
t18 = r_i_i_C(2) * t52;
t1 = [-t59 * t34 + t45 * t36 - t48 * t7 - t58 * t8, t39 * t36 + t19, t40 * t36 + t19 (-t32 * t51 + t34 * t35) * pkin(4) + t42, t42, t9; t58 * t10 + t45 * t34 + t59 * t36 + t48 * t9, t39 * t34 + t18, t40 * t34 + t18 (-t32 * t52 - t35 * t36) * pkin(4) + t43, t43, t7; 0, t29 + t41, t41 (-t46 - t57) * t26 + t47, -t26 * t46 + t47, t26 * t25;];
Ja_transl  = t1;
