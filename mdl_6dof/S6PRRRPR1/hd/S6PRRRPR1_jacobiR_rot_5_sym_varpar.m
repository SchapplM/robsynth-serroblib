% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:38
% EndTime: 2019-02-26 20:10:38
% DurationCPUTime: 0.04s
% Computational Cost: add. (89->14), mult. (87->32), div. (0->0), fcn. (134->8), ass. (0->26)
t98 = sin(pkin(11));
t99 = sin(pkin(6));
t109 = t98 * t99;
t100 = cos(pkin(11));
t108 = t100 * t99;
t103 = cos(qJ(2));
t107 = t103 * t99;
t102 = sin(qJ(2));
t106 = t99 * t102;
t101 = cos(pkin(6));
t105 = t101 * t102;
t104 = t101 * t103;
t97 = qJ(3) + qJ(4) + pkin(12);
t96 = cos(t97);
t95 = sin(t97);
t94 = t100 * t103 - t98 * t105;
t93 = -t100 * t102 - t98 * t104;
t92 = t100 * t105 + t98 * t103;
t91 = t100 * t104 - t98 * t102;
t90 = -t101 * t95 - t96 * t106;
t89 = t101 * t96 - t95 * t106;
t88 = -t95 * t109 - t94 * t96;
t87 = t96 * t109 - t94 * t95;
t86 = t95 * t108 - t92 * t96;
t85 = -t96 * t108 - t92 * t95;
t1 = [0, t93 * t96, t87, t87, 0, 0; 0, t91 * t96, t85, t85, 0, 0; 0, t96 * t107, t89, t89, 0, 0; 0, -t93 * t95, t88, t88, 0, 0; 0, -t91 * t95, t86, t86, 0, 0; 0, -t95 * t107, t90, t90, 0, 0; 0, t94, 0, 0, 0, 0; 0, t92, 0, 0, 0, 0; 0, t106, 0, 0, 0, 0;];
JR_rot  = t1;
