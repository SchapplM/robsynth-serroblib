% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR6_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:56:58
% EndTime: 2019-02-26 21:56:58
% DurationCPUTime: 0.04s
% Computational Cost: add. (68->18), mult. (76->8), div. (0->0), fcn. (122->6), ass. (0->14)
t52 = qJ(4) + qJ(5);
t50 = sin(t52);
t51 = cos(t52);
t53 = sin(qJ(2));
t55 = cos(qJ(2));
t60 = t55 * t50 - t53 * t51;
t57 = t53 * t50 + t55 * t51;
t56 = cos(qJ(1));
t54 = sin(qJ(1));
t44 = t57 * t56;
t43 = t60 * t56;
t42 = t57 * t54;
t41 = t60 * t54;
t1 = [-t42, t43, 0, -t43, -t43, 0; t44, t41, 0, -t41, -t41, 0; 0, t57, 0, -t57, -t57, 0; t41, t44, 0, -t44, -t44, 0; -t43, t42, 0, -t42, -t42, 0; 0, -t60, 0, t60, t60, 0; -t56, 0, 0, 0, 0, 0; -t54, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
