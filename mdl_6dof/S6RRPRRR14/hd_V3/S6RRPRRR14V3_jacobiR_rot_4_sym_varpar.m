% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR14V3_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobiR_rot_4_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:04
% EndTime: 2019-04-12 15:12:04
% DurationCPUTime: 0.03s
% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
t59 = sin(qJ(2));
t60 = sin(qJ(1));
t70 = t60 * t59;
t61 = cos(qJ(4));
t69 = t60 * t61;
t58 = sin(qJ(4));
t62 = cos(qJ(2));
t68 = t62 * t58;
t67 = t62 * t61;
t63 = cos(qJ(1));
t66 = t63 * t59;
t65 = t63 * t61;
t64 = t63 * t62;
t57 = t60 * t58 + t61 * t64;
t56 = -t58 * t64 + t69;
t55 = t63 * t58 - t60 * t67;
t54 = t60 * t68 + t65;
t1 = [t55, -t59 * t65, 0, t56, 0, 0; t57, -t59 * t69, 0, -t54, 0, 0; 0, t67, 0, -t59 * t58, 0, 0; t54, t58 * t66, 0, -t57, 0, 0; t56, t58 * t70, 0, t55, 0, 0; 0, -t68, 0, -t59 * t61, 0, 0; -t70, t64, 0, 0, 0, 0; t66, t60 * t62, 0, 0, 0, 0; 0, t59, 0, 0, 0, 0;];
JR_rot  = t1;
