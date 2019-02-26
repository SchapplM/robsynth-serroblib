% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPPR2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:39:45
% EndTime: 2019-02-26 20:39:45
% DurationCPUTime: 0.02s
% Computational Cost: add. (23->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->11)
t55 = qJ(3) + pkin(10);
t51 = sin(t55);
t56 = qJ(1) + pkin(9);
t52 = sin(t56);
t58 = t52 * t51;
t53 = cos(t55);
t54 = cos(t56);
t57 = t54 * t53;
t50 = t54 * t51;
t49 = t52 * t53;
t1 = [t54, 0, 0, 0, 0, 0; t52, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t49, 0, t50, 0, 0, 0; -t57, 0, t58, 0, 0, 0; 0, 0, -t53, 0, 0, 0; -t58, 0, t57, 0, 0, 0; t50, 0, t49, 0, 0, 0; 0, 0, t51, 0, 0, 0;];
JR_rot  = t1;
