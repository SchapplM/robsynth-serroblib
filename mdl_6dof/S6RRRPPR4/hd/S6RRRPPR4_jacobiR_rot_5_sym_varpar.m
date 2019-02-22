% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:51
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:51:32
% EndTime: 2019-02-22 11:51:32
% DurationCPUTime: 0.04s
% Computational Cost: add. (38->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
t85 = sin(qJ(2));
t86 = sin(qJ(1));
t93 = t86 * t85;
t84 = qJ(3) + pkin(10);
t82 = sin(t84);
t87 = cos(qJ(2));
t92 = t87 * t82;
t83 = cos(t84);
t91 = t87 * t83;
t88 = cos(qJ(1));
t90 = t88 * t85;
t89 = t88 * t87;
t81 = t86 * t82 + t83 * t89;
t80 = t82 * t89 - t86 * t83;
t79 = -t88 * t82 + t86 * t91;
t78 = -t88 * t83 - t86 * t92;
t1 = [-t79, -t83 * t90, -t80, 0, 0, 0; t81, -t83 * t93, t78, 0, 0, 0; 0, t91, -t85 * t82, 0, 0, 0; -t93, t89, 0, 0, 0, 0; t90, t86 * t87, 0, 0, 0, 0; 0, t85, 0, 0, 0, 0; t78, -t82 * t90, t81, 0, 0, 0; t80, -t82 * t93, t79, 0, 0, 0; 0, t92, t85 * t83, 0, 0, 0;];
JR_rot  = t1;
