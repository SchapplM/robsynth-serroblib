% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP1
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
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:39:33
% EndTime: 2019-02-26 22:39:33
% DurationCPUTime: 0.15s
% Computational Cost: add. (248->40), mult. (191->48), div. (0->0), fcn. (202->10), ass. (0->35)
t22 = qJ(2) + qJ(3);
t19 = qJ(4) + t22;
t15 = sin(t19);
t16 = cos(t19);
t24 = sin(qJ(5));
t44 = r_i_i_C(2) * t24;
t50 = r_i_i_C(3) * t16 + t15 * t44;
t26 = cos(qJ(5));
t17 = t26 * pkin(5) + pkin(4);
t23 = -qJ(6) - pkin(10);
t49 = t16 * t17 + (r_i_i_C(3) - t23) * t15;
t48 = pkin(5) + r_i_i_C(1);
t45 = r_i_i_C(1) * t26;
t32 = -t16 * t23 + (-t17 - t45) * t15;
t28 = t32 - pkin(3) * sin(t22);
t14 = pkin(3) * cos(t22);
t20 = cos(qJ(2)) * pkin(2);
t47 = pkin(1) + t14 + t20 + t49;
t25 = sin(qJ(1));
t41 = t50 * t25;
t27 = cos(qJ(1));
t40 = t50 * t27;
t39 = t25 * t24;
t38 = t25 * t26;
t37 = t27 * t24;
t36 = t27 * t26;
t34 = pkin(5) * t24 + pkin(7) + pkin(8) + pkin(9);
t3 = -t16 * t37 + t38;
t1 = t16 * t39 + t36;
t31 = (-t44 + t45) * t16 + t49;
t30 = -sin(qJ(2)) * pkin(2) + t28;
t29 = t14 + t31;
t4 = t16 * t36 + t39;
t2 = -t16 * t38 + t37;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t47 * t25 + t34 * t27, t30 * t27 + t40, t28 * t27 + t40, t32 * t27 + t40, -t4 * r_i_i_C(2) + t48 * t3, t27 * t15; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t34 * t25 + t47 * t27, t30 * t25 + t41, t28 * t25 + t41, t32 * t25 + t41, t2 * r_i_i_C(2) - t48 * t1, t25 * t15; 0, t20 + t29, t29, t31 (-r_i_i_C(2) * t26 - t48 * t24) * t15, -t16;];
Ja_transl  = t5;
