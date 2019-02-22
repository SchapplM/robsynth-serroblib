% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:09
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRP2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:09:47
% EndTime: 2019-02-22 11:09:47
% DurationCPUTime: 0.06s
% Computational Cost: add. (33->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
t79 = sin(qJ(5));
t80 = sin(qJ(1));
t86 = t80 * t79;
t81 = cos(qJ(5));
t85 = t80 * t81;
t78 = qJ(2) + pkin(9);
t76 = sin(t78);
t82 = cos(qJ(1));
t84 = t82 * t76;
t77 = cos(t78);
t83 = t82 * t77;
t75 = -t76 * t86 + t81 * t82;
t74 = t76 * t85 + t79 * t82;
t73 = t79 * t84 + t85;
t72 = t81 * t84 - t86;
t1 = [t75, t79 * t83, 0, 0, t72, 0; t73, t77 * t86, 0, 0, t74, 0; 0, t76 * t79, 0, 0, -t77 * t81, 0; -t74, t81 * t83, 0, 0, -t73, 0; t72, t77 * t85, 0, 0, t75, 0; 0, t76 * t81, 0, 0, t77 * t79, 0; -t80 * t77, -t84, 0, 0, 0, 0; t83, -t80 * t76, 0, 0, 0, 0; 0, t77, 0, 0, 0, 0;];
JR_rot  = t1;
