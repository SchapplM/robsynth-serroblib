% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR5_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobiR_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:48
% EndTime: 2019-02-26 21:30:48
% DurationCPUTime: 0.03s
% Computational Cost: add. (10->8), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
t55 = sin(qJ(2));
t56 = sin(qJ(1));
t62 = t56 * t55;
t57 = cos(qJ(2));
t61 = t56 * t57;
t58 = cos(qJ(1));
t60 = t58 * t55;
t59 = t58 * t57;
t54 = cos(pkin(6));
t53 = sin(pkin(6));
t52 = -t54 * t62 + t59;
t51 = t54 * t61 + t60;
t50 = t54 * t60 + t61;
t49 = t54 * t59 - t62;
t1 = [-t50, -t51, 0, 0, 0, 0; t52, t49, 0, 0, 0, 0; 0, t53 * t57, 0, 0, 0, 0; t49, t52, 0, 0, 0, 0; t51, t50, 0, 0, 0, 0; 0, t53 * t55, 0, 0, 0, 0; -t58 * t53, 0, 0, 0, 0, 0; -t56 * t53, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
