% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:46
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRP4_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_jacobiR_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:46:46
% EndTime: 2019-02-22 09:46:46
% DurationCPUTime: 0.08s
% Computational Cost: add. (19->13), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
t81 = sin(pkin(6));
t84 = sin(qJ(3));
t93 = t81 * t84;
t85 = sin(qJ(2));
t92 = t81 * t85;
t86 = cos(qJ(3));
t91 = t81 * t86;
t87 = cos(qJ(2));
t90 = t81 * t87;
t83 = cos(pkin(6));
t89 = t83 * t85;
t88 = t83 * t87;
t82 = cos(pkin(10));
t80 = sin(pkin(10));
t79 = -t80 * t89 + t82 * t87;
t78 = -t80 * t88 - t82 * t85;
t77 = t80 * t87 + t82 * t89;
t76 = -t80 * t85 + t82 * t88;
t1 = [0, t79, 0, 0, 0, 0; 0, t77, 0, 0, 0, 0; 0, t92, 0, 0, 0, 0; 0, -t78 * t86, t79 * t84 - t80 * t91, 0, 0, 0; 0, -t76 * t86, t77 * t84 + t82 * t91, 0, 0, 0; 0, -t86 * t90, -t83 * t86 + t84 * t92, 0, 0, 0; 0, t78 * t84, t79 * t86 + t80 * t93, 0, 0, 0; 0, t76 * t84, t77 * t86 - t82 * t93, 0, 0, 0; 0, t84 * t90, t83 * t84 + t85 * t91, 0, 0, 0;];
JR_rot  = t1;
