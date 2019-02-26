% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:44
% EndTime: 2019-02-26 20:08:44
% DurationCPUTime: 0.21s
% Computational Cost: add. (275->53), mult. (521->93), div. (0->0), fcn. (664->12), ass. (0->43)
t38 = cos(qJ(4));
t28 = t38 * pkin(4) + pkin(3);
t36 = sin(qJ(3));
t39 = cos(qJ(3));
t57 = r_i_i_C(2) + qJ(5) + pkin(9);
t59 = t28 * t39 + t57 * t36 + pkin(2);
t58 = pkin(5) + r_i_i_C(1);
t31 = qJ(4) + pkin(11);
t29 = sin(t31);
t56 = t29 * t39;
t30 = cos(t31);
t55 = t30 * t39;
t33 = sin(pkin(6));
t54 = t33 * t36;
t53 = t33 * t39;
t40 = cos(qJ(2));
t52 = t33 * t40;
t51 = t39 * t40;
t50 = r_i_i_C(3) + qJ(6);
t49 = cos(pkin(6));
t48 = cos(pkin(10));
t35 = sin(qJ(4));
t47 = pkin(4) * t35 + pkin(8);
t32 = sin(pkin(10));
t45 = t32 * t49;
t44 = t33 * t48;
t43 = t49 * t48;
t41 = -t50 * t29 - t58 * t30 - t28;
t37 = sin(qJ(2));
t24 = t49 * t36 + t37 * t53;
t23 = t37 * t54 - t49 * t39;
t22 = -t37 * t45 + t48 * t40;
t21 = t48 * t37 + t40 * t45;
t20 = t32 * t40 + t37 * t43;
t19 = t32 * t37 - t40 * t43;
t14 = t22 * t39 + t32 * t54;
t13 = t22 * t36 - t32 * t53;
t12 = t20 * t39 - t36 * t44;
t11 = t20 * t36 + t39 * t44;
t9 = t24 * t29 + t30 * t52;
t3 = t14 * t29 - t21 * t30;
t1 = t12 * t29 - t19 * t30;
t2 = [0, t58 * (-t21 * t55 + t22 * t29) + t50 * (-t21 * t56 - t22 * t30) + t47 * t22 - t59 * t21, t41 * t13 + t57 * t14, t50 * (t14 * t30 + t21 * t29) - t58 * t3 + (-t14 * t35 + t21 * t38) * pkin(4), t13, t3; 0, t58 * (-t19 * t55 + t20 * t29) + t50 * (-t19 * t56 - t20 * t30) + t47 * t20 - t59 * t19, t41 * t11 + t57 * t12, t50 * (t12 * t30 + t19 * t29) - t58 * t1 + (-t12 * t35 + t19 * t38) * pkin(4), t11, t1; 1 (t58 * (t29 * t37 + t30 * t51) + t50 * (t29 * t51 - t30 * t37) + t47 * t37 + t59 * t40) * t33, t41 * t23 + t57 * t24, -t58 * t9 + t50 * (t24 * t30 - t29 * t52) + (-t24 * t35 - t38 * t52) * pkin(4), t23, t9;];
Ja_transl  = t2;
