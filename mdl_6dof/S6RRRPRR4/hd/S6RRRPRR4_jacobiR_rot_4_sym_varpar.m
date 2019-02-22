% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:05
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR4_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_jacobiR_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:05:15
% EndTime: 2019-02-22 12:05:15
% DurationCPUTime: 0.04s
% Computational Cost: add. (35->12), mult. (38->18), div. (0->0), fcn. (66->6), ass. (0->20)
t87 = qJ(2) + qJ(3);
t86 = cos(t87);
t88 = sin(pkin(11));
t98 = t86 * t88;
t90 = sin(qJ(1));
t97 = t90 * t88;
t89 = cos(pkin(11));
t96 = t90 * t89;
t91 = cos(qJ(1));
t95 = t91 * t88;
t94 = t91 * t89;
t85 = sin(t87);
t93 = t85 * t96;
t92 = t85 * t94;
t84 = t91 * t86;
t83 = t90 * t86;
t82 = t86 * t89;
t81 = t85 * t95;
t80 = t85 * t97;
t1 = [-t86 * t96 + t95, -t92, -t92, 0, 0, 0; t86 * t94 + t97, -t93, -t93, 0, 0, 0; 0, t82, t82, 0, 0, 0; t86 * t97 + t94, t81, t81, 0, 0, 0; -t86 * t95 + t96, t80, t80, 0, 0, 0; 0, -t98, -t98, 0, 0, 0; -t90 * t85, t84, t84, 0, 0, 0; t91 * t85, t83, t83, 0, 0, 0; 0, t85, t85, 0, 0, 0;];
JR_rot  = t1;
