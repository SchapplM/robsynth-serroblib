% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPPRR2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:13
% EndTime: 2019-02-26 19:45:13
% DurationCPUTime: 0.05s
% Computational Cost: add. (41->13), mult. (122->36), div. (0->0), fcn. (178->10), ass. (0->23)
t95 = sin(pkin(6));
t99 = sin(qJ(5));
t106 = t95 * t99;
t101 = cos(qJ(5));
t105 = t101 * t95;
t100 = sin(qJ(2));
t102 = cos(qJ(2));
t93 = sin(pkin(11));
t96 = cos(pkin(11));
t104 = t100 * t96 + t102 * t93;
t91 = t100 * t93 - t102 * t96;
t98 = cos(pkin(6));
t103 = t91 * t98;
t97 = cos(pkin(10));
t94 = sin(pkin(10));
t90 = t104 * t98;
t89 = t104 * t95;
t88 = t91 * t95;
t86 = t94 * t103 - t104 * t97;
t85 = -t94 * t90 - t97 * t91;
t84 = -t97 * t103 - t104 * t94;
t83 = t97 * t90 - t94 * t91;
t1 = [0, t85 * t99, 0, 0, -t86 * t101 - t94 * t106, 0; 0, t83 * t99, 0, 0, -t84 * t101 + t97 * t106, 0; 0, t89 * t99, 0, 0, t88 * t101 - t98 * t99, 0; 0, t85 * t101, 0, 0, -t94 * t105 + t86 * t99, 0; 0, t83 * t101, 0, 0, t97 * t105 + t84 * t99, 0; 0, t89 * t101, 0, 0, -t98 * t101 - t88 * t99, 0; 0, t86, 0, 0, 0, 0; 0, t84, 0, 0, 0, 0; 0, -t88, 0, 0, 0, 0;];
JR_rot  = t1;
