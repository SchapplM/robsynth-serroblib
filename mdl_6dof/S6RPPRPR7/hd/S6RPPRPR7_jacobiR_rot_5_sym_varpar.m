% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRPR7_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:29:09
% EndTime: 2019-02-26 20:29:09
% DurationCPUTime: 0.03s
% Computational Cost: add. (25->11), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
t54 = sin(pkin(10));
t56 = sin(qJ(1));
t61 = t56 * t54;
t55 = cos(pkin(10));
t60 = t56 * t55;
t57 = cos(qJ(1));
t59 = t57 * t54;
t58 = t57 * t55;
t53 = pkin(9) + qJ(4);
t52 = cos(t53);
t51 = sin(t53);
t1 = [t51 * t58 - t61, 0, 0, t52 * t60, 0, 0; t51 * t60 + t59, 0, 0, -t52 * t58, 0, 0; 0, 0, 0, -t51 * t55, 0, 0; -t51 * t59 - t60, 0, 0, -t52 * t61, 0, 0; -t51 * t61 + t58, 0, 0, t52 * t59, 0, 0; 0, 0, 0, t51 * t54, 0, 0; -t57 * t52, 0, 0, t56 * t51, 0, 0; -t56 * t52, 0, 0, -t57 * t51, 0, 0; 0, 0, 0, t52, 0, 0;];
JR_rot  = t1;
