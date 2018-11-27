% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JR_rot [9x7]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JR_rot = S7RRRRRRR1_jacobiR_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiR_rot_3_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiR_rot_3_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 21:20:53
% EndTime: 2018-11-26 21:20:53
% DurationCPUTime: 0.04s
% Computational Cost: add. (17->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
t58 = sin(qJ(2));
t59 = sin(qJ(1));
t69 = t59 * t58;
t60 = cos(qJ(3));
t68 = t59 * t60;
t57 = sin(qJ(3));
t61 = cos(qJ(2));
t67 = t61 * t57;
t66 = t61 * t60;
t62 = cos(qJ(1));
t65 = t62 * t58;
t64 = t62 * t60;
t63 = t62 * t61;
t56 = -t59 * t57 + t60 * t63;
t55 = -t57 * t63 - t68;
t54 = -t62 * t57 - t59 * t66;
t53 = t59 * t67 - t64;
t1 = [t54, -t58 * t64, t55, 0, 0, 0, 0; t56, -t58 * t68, -t53, 0, 0, 0, 0; 0, t66, -t58 * t57, 0, 0, 0, 0; t53, t57 * t65, -t56, 0, 0, 0, 0; t55, t57 * t69, t54, 0, 0, 0, 0; 0, -t67, -t58 * t60, 0, 0, 0, 0; t69, -t63, 0, 0, 0, 0, 0; -t65, -t59 * t61, 0, 0, 0, 0, 0; 0, -t58, 0, 0, 0, 0, 0;];
JR_rot  = t1;
