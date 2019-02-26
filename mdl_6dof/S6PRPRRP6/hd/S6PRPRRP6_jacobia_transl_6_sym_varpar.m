% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRP6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:15
% EndTime: 2019-02-26 19:53:15
% DurationCPUTime: 0.15s
% Computational Cost: add. (182->46), mult. (470->82), div. (0->0), fcn. (603->10), ass. (0->36)
t52 = pkin(2) + pkin(8);
t51 = pkin(5) + r_i_i_C(1);
t50 = pkin(9) + r_i_i_C(2);
t29 = sin(pkin(6));
t33 = sin(qJ(4));
t49 = t29 * t33;
t36 = cos(qJ(4));
t48 = t29 * t36;
t37 = cos(qJ(2));
t47 = t29 * t37;
t31 = cos(pkin(6));
t34 = sin(qJ(2));
t46 = t31 * t34;
t45 = t31 * t37;
t32 = sin(qJ(5));
t44 = t32 * t33;
t43 = t32 * t34;
t35 = cos(qJ(5));
t42 = t33 * t35;
t41 = t34 * t35;
t40 = r_i_i_C(3) + qJ(6);
t39 = pkin(4) * t33 - t50 * t36 + qJ(3);
t38 = t40 * t32 + t51 * t35 + pkin(4);
t30 = cos(pkin(10));
t28 = sin(pkin(10));
t25 = t31 * t36 - t33 * t47;
t23 = -t28 * t46 + t30 * t37;
t22 = t28 * t45 + t30 * t34;
t21 = t28 * t37 + t30 * t46;
t20 = t28 * t34 - t30 * t45;
t13 = t25 * t32 - t29 * t41;
t12 = -t20 * t33 + t30 * t48;
t10 = t22 * t33 + t28 * t48;
t3 = -t12 * t32 - t21 * t35;
t1 = t10 * t32 - t23 * t35;
t2 = [0, t51 * (-t22 * t32 + t23 * t42) + t40 * (t22 * t35 + t23 * t44) - t52 * t22 + t39 * t23, t22, t50 * t10 + t38 * (t22 * t36 - t28 * t49) t40 * (t10 * t35 + t23 * t32) - t51 * t1, t1; 0, t51 * (-t20 * t32 + t21 * t42) + t40 * (t20 * t35 + t21 * t44) - t52 * t20 + t39 * t21, t20, -t50 * t12 + t38 * (t20 * t36 + t30 * t49) t40 * (-t12 * t35 + t21 * t32) - t51 * t3, t3; 1 (t51 * (t32 * t37 + t33 * t41) + t40 * (t33 * t43 - t35 * t37) + t39 * t34 + t52 * t37) * t29, -t47, t50 * t25 + t38 * (-t31 * t33 - t36 * t47) t40 * (t25 * t35 + t29 * t43) - t51 * t13, t13;];
Ja_transl  = t2;
