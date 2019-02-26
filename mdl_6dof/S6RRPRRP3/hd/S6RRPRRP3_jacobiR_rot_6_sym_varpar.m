% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:47:10
% EndTime: 2019-02-26 21:47:10
% DurationCPUTime: 0.04s
% Computational Cost: add. (81->18), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
t96 = qJ(2) + pkin(10);
t92 = sin(t96);
t97 = qJ(4) + qJ(5);
t94 = sin(t97);
t105 = t92 * t94;
t95 = cos(t97);
t104 = t92 * t95;
t98 = sin(qJ(1));
t103 = t98 * t94;
t102 = t98 * t95;
t99 = cos(qJ(1));
t101 = t99 * t94;
t100 = t99 * t95;
t93 = cos(t96);
t91 = t93 * t100 + t103;
t90 = -t93 * t101 + t102;
t89 = -t93 * t102 + t101;
t88 = t93 * t103 + t100;
t1 = [t89, -t92 * t100, 0, t90, t90, 0; t91, -t92 * t102, 0, -t88, -t88, 0; 0, t93 * t95, 0, -t105, -t105, 0; t88, t92 * t101, 0, -t91, -t91, 0; t90, t92 * t103, 0, t89, t89, 0; 0, -t93 * t94, 0, -t104, -t104, 0; -t98 * t92, t99 * t93, 0, 0, 0, 0; t99 * t92, t98 * t93, 0, 0, 0, 0; 0, t92, 0, 0, 0, 0;];
JR_rot  = t1;
